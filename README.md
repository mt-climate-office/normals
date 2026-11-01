# normals 2.0 (in development)

`normals` 2.0 is a rewrite of this package. It holds the climate-normals methods used by the
Montana Climate Office's SNOTEL and SNODAS pipelines: empirical, fitted and smoothed normals,
distribution families, percentiles and standardized anomalies, the SWE outlook, and shared plot
helpers. It reads and writes no storage and has no AWS dependency.

The API design is in [docs/API.md](docs/API.md). The functions are stubs until it is approved.

## normals 1.x

The 1.x package (gamma normals from gridMET) is unchanged at the tags
[`v1.0.0`](https://github.com/mt-climate-office/normals/tree/v1.0.0) and
[`v1.0.1`](https://github.com/mt-climate-office/normals/tree/v1.0.1), and at commit `562af10`,
which is where this branch starts. Users of 1.x should pin one of those. If you use 1.x in a
publication, please cite:

    Colin, Brust, R. Kyle, Bocinsky, Zachary, Hoylman, & Kelsey, Jencso. (2022). Normals: An R
    Package for Creating Montana's Climate Normals (1.0.1). Zenodo.
    https://doi.org/10.5281/zenodo.7278462
