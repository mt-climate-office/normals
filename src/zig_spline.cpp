// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

// Smooth zero-inflated gamma (ZIG) normals: the daily ZIG (src/zig.cpp)
// smoothed over the year. Each day of the water year contributes the daily
// ZIG's own statistics (zig_lmom()): n (years observed), z (zero years), and
// the L-moments l1 and L-CV of the positive values. Three cyclic P-splines in
// day of water year are fit to them, per pixel:
//
//   logit(p_zero) = B beta_nu    binomial penalised likelihood on z of n
//                                (the smoothed form of the daily z / n)
//   log(l1)       = B beta_l1    weighted penalised least squares on the
//   log(L-CV)     = B beta_lcv   daily L-moments, weights = positive years
//
// Gamma shape follows from the smoothed L-CV by Hosking's approximation and
// scale = l1 / shape, exactly as zig_fit() does for one day
// (zig_spline_params() in R/calc_snodas_zig_spline_normals.R).
//
// B is a cyclic cubic B-spline basis (mgcv bs = "cp", no constraint
// absorbed): 4 nonzeros per row and B %*% rep(c, k) == c. The penalty is a
// cyclic difference penalty whose null space is the constants. Smoothing
// parameters are fixed.
//
// The L-moment curves are linear smoothers, solved directly. Only the
// binomial curve iterates (penalised Fisher scoring), from a start that is
// one weighted least-squares solve on the empirical logit.
//
// Moved from mco-snodas src/zig_spline.cpp (62fe19a); only the export's name
// changed (cpp_zig_spline_fit, internal). family_zig()'s smoother calls it.

namespace {

struct Basis {
  int n, k, m;          // rows (days), columns, max nonzeros per row
  arma::umat idx;       // n x m column index of each nonzero
  arma::mat  val;       // n x m value of each nonzero
  arma::mat  S;         // k x k penalty
};

Basis make_basis(const NumericMatrix& X, const NumericMatrix& S) {
  Basis b;
  b.n = X.nrow();
  b.k = X.ncol();
  b.m = 0;
  for (int i = 0; i < b.n; i++) {
    int c = 0;
    for (int j = 0; j < b.k; j++) if (X(i, j) != 0.0) c++;
    b.m = std::max(b.m, c);
  }
  b.idx.zeros(b.n, b.m);
  b.val.zeros(b.n, b.m);
  for (int i = 0; i < b.n; i++) {
    int c = 0;
    for (int j = 0; j < b.k; j++) {
      if (X(i, j) != 0.0) {
        b.idx(i, c) = j;
        b.val(i, c) = X(i, j);
        c++;
      }
    }
  }
  b.S = arma::symmatu(as<arma::mat>(S));
  return b;
}

inline void linear_predictor(const Basis& b, const arma::vec& beta, arma::vec& eta) {
  for (int i = 0; i < b.n; i++) {
    double e = 0.0;
    for (int c = 0; c < b.m; c++) e += b.val(i, c) * beta(b.idx(i, c));
    eta(i) = e;
  }
}

// H = X' diag(w) X, g = X' u, exploiting the sparse rows of X.
void accumulate(const Basis& b, const arma::vec& w, const arma::vec& u,
                arma::mat& H, arma::vec& g) {
  H.zeros();
  g.zeros();
  for (int i = 0; i < b.n; i++) {
    if (w(i) == 0.0 && u(i) == 0.0) continue;
    for (int c1 = 0; c1 < b.m; c1++) {
      const arma::uword j1 = b.idx(i, c1);
      const double v1 = b.val(i, c1);
      g(j1) += v1 * u(i);
      const double wv1 = w(i) * v1;
      for (int c2 = 0; c2 < b.m; c2++) H(j1, b.idx(i, c2)) += wv1 * b.val(i, c2);
    }
  }
  H = arma::symmatu(H);   // exact symmetry for the sympd solvers
}

inline double quad(const arma::mat& S, const arma::vec& beta) {
  return arma::as_scalar(beta.t() * S * beta);
}

inline double log1pexp(double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

const double RIDGE = 1e-8;   // keeps H + lambda S positive definite

// Solve (H + lambda S) beta = g.
bool penalised_solve(const Basis& b, const arma::mat& H, const arma::vec& g,
                     double lambda, arma::vec& beta) {
  arma::mat A = arma::symmatu(H + lambda * b.S);
  A.diag() += RIDGE;
  return arma::solve(beta, A, g, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);
}

// Effective degrees of freedom tr((H + lambda S)^-1 H).
double edf(const Basis& b, const arma::mat& H, double lambda) {
  arma::mat A = arma::symmatu(H + lambda * b.S);
  A.diag() += RIDGE;
  arma::mat M;
  if (!arma::solve(M, A, H, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx)) return NA_REAL;
  return arma::trace(M);
}

// Per-pixel daily statistics, one entry per day of the water year.
struct Day {
  arma::vec n, z;         // years observed, zero years
  arma::vec w;            // L-moment weight: positive years, 0 without L-moments
  arma::vec y_l1, y_lcv;  // log l1, log L-CV (0 where w = 0)
};

// ── nu: binomial on zero indicators ─────────────────────────────────────────
// Days with no snow in any year are quasi-separated (MLE logit(nu) -> +inf),
// so logit(nu) is softly bounded to +/- ETA_NU_MAX by a quadratic barrier,
// -0.5 K (eta - bound)^2 outside the bounds: p_zero within 5e-5 of 0 or 1,
// finer than ~25 years can resolve. Unlike a hard wall the barrier never
// stalls step halving at the bound.

const double ETA_NU_MAX = 10.0;
const double BARRIER_K  = 100.0;

inline double barrier(double e) {
  const double x = e > ETA_NU_MAX ? e - ETA_NU_MAX : (e < -ETA_NU_MAX ? e + ETA_NU_MAX : 0.0);
  return -0.5 * BARRIER_K * x * x;
}

inline void barrier_score(double e, double& u, double& w) {
  if (e > ETA_NU_MAX)       { u -= BARRIER_K * (e - ETA_NU_MAX); w += BARRIER_K; }
  else if (e < -ETA_NU_MAX) { u -= BARRIER_K * (e + ETA_NU_MAX); w += BARRIER_K; }
}

double ll_nu(const Basis& b, const Day& d, const arma::vec& beta, arma::vec& eta) {
  linear_predictor(b, beta, eta);
  double ll = 0.0;
  for (int i = 0; i < b.n; i++) {
    if (d.n(i) > 0.0) ll += d.z(i) * eta(i) - d.n(i) * log1pexp(eta(i));
    ll += barrier(eta(i));
  }
  return ll;
}

// with_barrier = false gives the likelihood-only weights used for the EDF.
void score_nu(const Basis& b, const Day& d, const arma::vec& beta, arma::vec& eta,
              arma::vec& w, arma::vec& u, bool with_barrier = true) {
  linear_predictor(b, beta, eta);
  for (int i = 0; i < b.n; i++) {
    const double pr = 1.0 / (1.0 + std::exp(-eta(i)));
    u(i) = d.z(i) - d.n(i) * pr;
    w(i) = d.n(i) * pr * (1.0 - pr);
    if (with_barrier) barrier_score(eta(i), u(i), w(i));
  }
}

// Start: weighted penalised least squares on the empirical logit,
// log((z + 0.5) / (n - z + 0.5)), with binomial weights. This is one scoring
// step from the saturated fit; the 0.5 only places the start, the optimum is
// the binomial one.
void start_nu(const Basis& b, const Day& d, double lambda, arma::vec& beta) {
  arma::vec w(b.n), u(b.n), g(b.k);
  arma::mat H(b.k, b.k);
  for (int i = 0; i < b.n; i++) {
    if (d.n(i) <= 0.0) { w(i) = 0.0; u(i) = 0.0; continue; }
    const double p = (d.z(i) + 0.5) / (d.n(i) + 1.0);
    const double e = std::min(ETA_NU_MAX, std::max(-ETA_NU_MAX, std::log(p / (1.0 - p))));
    w(i) = d.n(i) * p * (1.0 - p);
    u(i) = w(i) * e;
  }
  accumulate(b, w, u, H, g);
  if (!penalised_solve(b, H, g, lambda, beta)) beta.zeros();
}

struct FitOut {
  bool   converged;
  int    iter;
  double edf;
};

// Penalised Fisher scoring with step halving; converged when the Newton
// decrement (twice the predicted gain of a full step) falls below tol, or
// when no step improves the penalised likelihood (optimal to floating-point
// precision).
FitOut fit_nu(const Basis& b, const Day& d, arma::vec& beta, double lambda,
              int maxit, double tol) {
  arma::vec eta(b.n), w(b.n), u(b.n), g(b.k), delta;
  arma::mat H(b.k, b.k);
  auto pll = [&](const arma::vec& be) {
    arma::vec e(b.n);
    return ll_nu(b, d, be, e) - 0.5 * lambda * quad(b.S, be);
  };
  FitOut out{false, 0, NA_REAL};
  for (int it = 1; it <= maxit; it++) {
    out.iter = it;
    score_nu(b, d, beta, eta, w, u);
    accumulate(b, w, u, H, g);
    const arma::vec grad = g - lambda * (b.S * beta);
    if (!penalised_solve(b, H, grad, lambda, delta)) break;
    const double decrement = arma::dot(grad, delta);
    const double pll0 = pll(beta);
    bool moved = false;
    double step = 1.0;
    for (int h = 0; h < 30; h++) {
      const arma::vec cand = beta + step * delta;
      const double pl = pll(cand);
      if (std::isfinite(pl) && pl >= pll0) { beta = cand; moved = true; break; }
      step *= 0.5;
    }
    if (!moved || decrement < tol) { out.converged = true; break; }
  }
  score_nu(b, d, beta, eta, w, u, false);
  accumulate(b, w, u, H, g);
  out.edf = edf(b, H, lambda);
  return out;
}

// ── l1, L-CV: weighted penalised least squares on the log scale ─────────────
// One X'WX per pixel serves both curves. Days without L-moments (fewer than
// 4 positive years) have weight 0, and the curves bridge them.

void fit_lmom(const Basis& b, const Day& d, double lambda_l1, double lambda_lcv,
              arma::vec& beta_l1, arma::vec& beta_lcv, double& edf_l1, double& edf_lcv) {
  arma::vec g1(b.k), g2(b.k), wy(b.n);
  arma::mat H(b.k, b.k);
  for (int i = 0; i < b.n; i++) wy(i) = d.w(i) * d.y_l1(i);
  accumulate(b, d.w, wy, H, g1);
  for (int i = 0; i < b.n; i++) wy(i) = d.w(i) * d.y_lcv(i);
  accumulate(b, d.w, wy, H, g2);
  if (!penalised_solve(b, H, g1, lambda_l1, beta_l1))   beta_l1.fill(NA_REAL);
  if (!penalised_solve(b, H, g2, lambda_lcv, beta_lcv)) beta_lcv.fill(NA_REAL);
  edf_l1  = edf(b, H, lambda_l1);
  edf_lcv = edf(b, H, lambda_lcv);
}

} // namespace


// [[Rcpp::export]]
List cpp_zig_spline_fit(NumericMatrix n_obs,
                        NumericMatrix n_zero,
                        NumericMatrix l1,
                        NumericMatrix lcv,
                        NumericMatrix basis,
                        NumericMatrix penalty,
                        NumericVector lambda,
                        int    min_pos  = 10,
                        int    min_days = 30,
                        int    maxit    = 100,
                        double tol      = 1e-6) {
  // n_obs, n_zero, l1, lcv: n_pix x n_day daily statistics (zig_lmom() per
  //   day); a row with any NA n_obs is status 3. l1 and lcv are NA on days
  //   without L-moments.
  // basis:    n_day x k cyclic B-spline basis; penalty: k x k.
  // lambda:   smoothing parameters (nu, l1, lcv).
  // min_pos:  fewer positive values on record → fit nu only.
  // min_days: fewer days with L-moments → l1 and L-CV constant over the year
  //           (their weighted means on the log scale); none → nu only.
  // tol:      the nu fit has converged when the Newton decrement falls below this.
  //
  // Returns a named list:
  //   $coef   n_pix x 3k  [beta_nu | beta_l1 | beta_lcv]
  //   $edf    n_pix x 3   effective degrees of freedom
  //   $iter   n_pix       scoring iterations of the nu fit
  //   $status n_pix       0 ok · 1 no snow on record · 2 nu only (too few
  //                       positives or no day with L-moments) · 3 no data ·
  //                       4 nu not converged (l1 and L-CV still fit) ·
  //                       5 ok, l1 and L-CV constant

  const int n_pix = n_obs.nrow();
  const int n_day = n_obs.ncol();
  if (basis.nrow() != n_day) stop("basis must have one row per day");
  if (lambda.size() != 3) stop("lambda must be (nu, l1, lcv)");

  const Basis b = make_basis(basis, penalty);
  const int k = b.k;

  NumericMatrix coef_out(n_pix, 3 * k), edf_out(n_pix, 3);
  NumericVector iter_out(n_pix, NA_REAL);
  IntegerVector status(n_pix);
  std::fill(coef_out.begin(), coef_out.end(), NA_REAL);
  std::fill(edf_out.begin(), edf_out.end(), NA_REAL);

  Day d;
  d.n.set_size(n_day); d.z.set_size(n_day); d.w.set_size(n_day);
  d.y_l1.set_size(n_day); d.y_lcv.set_size(n_day);
  arma::vec beta_nu(k), beta_l1(k), beta_lcv(k);

  for (int i = 0; i < n_pix; i++) {
    if ((i + 1) % 10000 == 0) checkUserInterrupt();

    bool any_na = false;
    int n_lmom = 0;
    for (int j = 0; j < n_day; j++) {
      const double nn = n_obs(i, j);
      if (ISNAN(nn)) { any_na = true; break; }
      d.n(j) = nn;
      d.z(j) = n_zero(i, j);
      const bool ok = !ISNAN(l1(i, j)) && !ISNAN(lcv(i, j)) && l1(i, j) > 0.0 && lcv(i, j) > 0.0;
      d.w(j)     = ok ? nn - n_zero(i, j) : 0.0;
      d.y_l1(j)  = ok ? std::log(l1(i, j))  : 0.0;
      d.y_lcv(j) = ok ? std::log(lcv(i, j)) : 0.0;
      if (ok) n_lmom++;
    }
    const double N = any_na ? 0.0 : arma::accu(d.n);
    const double P = any_na ? 0.0 : N - arma::accu(d.z);
    if (N == 0.0) { status[i] = 3; continue; }
    if (P == 0.0) { status[i] = 1; continue; }

    // ── nu ──
    start_nu(b, d, lambda[0], beta_nu);
    const FitOut fn = fit_nu(b, d, beta_nu, lambda[0], maxit, tol);
    for (int j = 0; j < k; j++) coef_out(i, j) = beta_nu(j);
    edf_out(i, 0) = fn.edf;
    iter_out[i]   = fn.iter;
    if (P < min_pos || n_lmom == 0) { status[i] = fn.converged ? 2 : 4; continue; }

    // ── l1, L-CV ──
    if (n_lmom < min_days) {
      // Constant curves: B-splines sum to one, so beta = c * 1 gives eta = c.
      const double W = arma::accu(d.w);
      beta_l1.fill(arma::dot(d.w, d.y_l1) / W);
      beta_lcv.fill(arma::dot(d.w, d.y_lcv) / W);
      edf_out(i, 1) = edf_out(i, 2) = 1.0;
    } else {
      fit_lmom(b, d, lambda[1], lambda[2], beta_l1, beta_lcv, edf_out(i, 1), edf_out(i, 2));
    }
    status[i] = !fn.converged ? 4 : n_lmom < min_days ? 5 : 0;
    for (int j = 0; j < k; j++) {
      coef_out(i, k + j)     = beta_l1(j);
      coef_out(i, 2 * k + j) = beta_lcv(j);
    }
  }

  return List::create(
    Named("coef")   = coef_out,
    Named("edf")    = edf_out,
    Named("iter")   = iter_out,
    Named("status") = status
  );
}
