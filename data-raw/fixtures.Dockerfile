# The pinned fixture platform (API.md section 7): the image the golden
# fixtures are made in, and whose base CI's R-CMD-check job runs in.
#
#   docker build --platform linux/amd64 -f data-raw/fixtures.Dockerfile -t normals-fixtures data-raw
#
# Base: rocker/r-ver 4.5.3 (the mco-snotel base image), by digest. CRAN
# packages come from its P3M snapshot (2026-04-23, set in the image). BLAS
# and LAPACK are the image's OpenBLAS 0.3.26, single-threaded and with its
# kernel pinned (OPENBLAS_CORETYPE), so every x86-64 machine computes the
# same bits.
FROM rocker/r-ver:4.5.3@sha256:f5a72ffeea7ed592b2a5990f87101beef8ef3cae0e4ff7327088e80368f0aef9

ENV OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_CORETYPE=Nehalem

# Runtime libraries of the P3M binaries below, and git (provenance checks).
RUN apt-get update && apt-get install -y --no-install-recommends \
      git libgdal34t64 libgeos-c1t64 libproj25 libudunits2-0 \
      libpng16-16t64 libfreetype6 libfontconfig1 libharfbuzz0b libfribidi0 \
      libtiff6 libjpeg-turbo8 libwebp7 libwebpmux3 libcurl4t64 libssl3t64 \
  && rm -rf /var/lib/apt/lists/*

# The package's dependencies and Suggests, and what the pipelines' code
# needs to make the fixtures (arrow, terra, dplyr, ...).
RUN Rscript -e 'install.packages(c( \
      "Rcpp", "RcppArmadillo", "tibble", "testthat", "lubridate", \
      "fastpng", "ggplot2", "gtable", "patchwork", "png", "ragg", "scales", "systemfonts", \
      "arrow", "terra", "units", "dplyr", "dtplyr", "purrr", "magrittr", "jsonlite", "tidyr", "stringr"))' \
  && Rscript -e 'for (p in c("Rcpp", "RcppArmadillo", "arrow", "terra", "units", "ragg", "systemfonts", "testthat")) library(p, character.only = TRUE)'
