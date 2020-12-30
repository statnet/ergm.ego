# `ergm.ego`: Fit, Simulate and Diagnose Exponential-Family Random Graph Models to Egocentrically Sampled Network Data

[![Build Status](https://travis-ci.org/statnet/ergm.ego.svg?branch=master)](https://travis-ci.org/statnet/ergm.ego)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ergm.ego?color=2ED968)](https://cranlogs.r-pkg.org/)
[![cran version](https://www.r-pkg.org/badges/version/ergm.ego)](https://cran.r-project.org/package=ergm.ego)
[![Coverage status](https://codecov.io/gh/statnet/ergm.ego/branch/master/graph/badge.svg)](https://codecov.io/github/statnet/ergm.ego?branch=master)
[![R build status](https://github.com/statnet/ergm.ego/workflows/R-CMD-check/badge.svg)](https://github.com/statnet/ergm.ego/actions)

Utilities for managing egocentrically sampled network data and a wrapper around the 'ergm' package to facilitate ERGM inference and simulation from such data.

## Public and Private repositories

To facilitate open development of the package while giving the core developers an opportunity to publish on their developments before opening them up for general use, this project comprises two repositories:
* A public repository `statnet/ergm.ego`
* A private repository `statnet/ergm.ego-private`

The intention is that all developments in `statnet/ergm.ego-private` will eventually make their way into `statnet/ergm.ego` and onto CRAN.

Developers and Contributing Users to the Statnet Project should read https://statnet.github.io/private/ for information about the relationship between the public and the private repository and the workflows involved.

## Latest Windows and MacOS binaries

A set of binaries is built after every commit to the repository. We strongly encourage testing against them before filing a bug report, as they may contain fixes that have not yet been sent to CRAN. They can be downloaded through the following links:

* [MacOS binary (a `.tgz` file in a `.zip` file)](https://nightly.link/statnet/ergm.ego/workflows/R-CMD-check.yaml/master/macOS-rrelease-binaries.zip)
* [Windows binary (a `.zip` file in a `.zip` file)](https://nightly.link/statnet/ergm.ego/workflows/R-CMD-check.yaml/master/Windows-rrelease-binaries.zip)

You will need to extract the MacOS `.tgz` or the Windows `.zip` file from the outer `.zip` file before installing. These binaries are usually built under the latest version of R and their operating system and may not work under other versions.

You may also want to install the corresponding latest binaries for packages on which `ergm.ego` depends, in particular [`statnet.common`](https://github.com/statnet/statnet.common), [`rle`](https://github.com/statnet/rle), [`network`](https://github.com/statnet/network), and [`ergm`](https://github.com/statnet/ergm).
