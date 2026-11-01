#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Standard normal quantiles for grids of percentiles: percentile p in [0, 100]
// -> qnorm(p / 100) clamped to [-4, 4], x 100, rounded to an integer (e.g.
// -196 for z = -1.96; p <= 0 -> -400, p >= 100 -> 400, NA -> NA). Stored as
// INT2S; divide by 100 for z. Uses Acklam's rational approximation (relative
// error ~1.15e-9).
//
// Moved from mco-snodas src/qnorm_grid.cpp (62fe19a). Changes: the export is
// internal (cpp_qnorm_grid) and the clamp is the `z_max` argument (4 today).

// Acklam's coefficients.
static const double A1 = -3.969683028665376e+01;
static const double A2 =  2.209460984245205e+02;
static const double A3 = -2.759285104469687e+02;
static const double A4 =  1.383577518672690e+02;
static const double A5 = -3.066479806614716e+01;
static const double A6 =  2.506628277459239e+00;

static const double B1 = -5.447609879822406e+01;
static const double B2 =  1.615858368580409e+02;
static const double B3 = -1.556989798598866e+02;
static const double B4 =  6.680131188771972e+01;
static const double B5 = -1.328068155288572e+01;

static const double C1 = -7.784894002430293e-03;
static const double C2 = -3.223964580411365e-01;
static const double C3 = -2.400758277161838e+00;
static const double C4 = -2.549732539343734e+00;
static const double C5 =  4.374664141464968e+00;
static const double C6 =  2.938163982698783e+00;

static const double D1 =  7.784695709041462e-03;
static const double D2 =  3.224671290700398e-01;
static const double D3 =  2.445134137142996e+00;
static const double D4 =  3.754408661907416e+00;

static const double P_LOW  = 0.02425;
static const double P_HIGH = 1.0 - P_LOW;

inline double qnorm_scalar(double p) {
  // p on the [0, 1] scale.
  if (p <= 0.0) return R_NegInf;
  if (p >= 1.0) return R_PosInf;

  double q, r;

  if (p < P_LOW) {
    // Lower tail
    q = std::sqrt(-2.0 * std::log(p));
    return (((((C1*q + C2)*q + C3)*q + C4)*q + C5)*q + C6) /
            ((((D1*q + D2)*q + D3)*q + D4)*q + 1.0);
  } else if (p <= P_HIGH) {
    // Central region
    q = p - 0.5;
    r = q * q;
    return (((((A1*r + A2)*r + A3)*r + A4)*r + A5)*r + A6) * q /
           (((((B1*r + B2)*r + B3)*r + B4)*r + B5)*r + 1.0);
  } else {
    // Upper tail
    q = std::sqrt(-2.0 * std::log(1.0 - p));
    return -(((((C1*q + C2)*q + C3)*q + C4)*q + C5)*q + C6) /
             ((((D1*q + D2)*q + D3)*q + D4)*q + 1.0);
  }
}


// [[Rcpp::export]]
IntegerVector cpp_qnorm_grid(NumericVector ptile, double z_max) {
  // ptile: percentiles on the 0-100 scale; z_max: the clamp.

  int n = ptile.size();
  IntegerVector out(n);

  for (int i = 0; i < n; i++) {
    double v = ptile[i];
    if (ISNAN(v)) {
      out[i] = NA_INTEGER;
    } else {
      double z = qnorm_scalar(v / 100.0);
      z = std::max(-z_max, std::min(z_max, z));
      out[i] = (int)std::round(z * 100.0);
    }
  }

  return out;
}
