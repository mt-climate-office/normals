# The fixture platform's fingerprint: R, the OS architecture, the BLAS and
# LAPACK R links, and the OpenBLAS kernel and thread settings. The golden
# fixtures are made on one pinned platform (data-raw/fixtures.Dockerfile);
# results that pass through linear algebra are compared only where this
# fingerprint matches the fixtures' (tests/testthat/fixtures/platform.rds).
# data-raw/make-fixtures.R sources this file.
golden_fingerprint <- function() {
  si <- utils::sessionInfo()
  list(R = R.version.string, platform = R.version$platform,
       blas = si$BLAS, lapack = si$LAPACK,
       openblas_coretype = Sys.getenv("OPENBLAS_CORETYPE"),
       openblas_num_threads = Sys.getenv("OPENBLAS_NUM_THREADS"),
       omp_num_threads = Sys.getenv("OMP_NUM_THREADS"))
}

# A file's SHA-256 (tools::sha256sum() from R 4.5).
golden_sha256 <- function(path) {
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE))
    return(unname(tools::sha256sum(path)))
  sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
}
