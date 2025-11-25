# spatstat.random

## Random Generation and Simulation for the spatstat family

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.random)](http://CRAN.R-project.org/package=spatstat.random) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.random)](https://github.com/spatstat/spatstat.random)

The original `spatstat` package has been split into
several sub-packages
(see [spatstat/spatstat](https://github.com/spatstat/spatstat))

This package `spatstat.random` is one of the sub-packages.
It contains the functions for **random generation** of data
and **simulation** of models.

You are viewing the GitHub repository which holds
the latest **development version** of `spatstat.random`.
For the latest public release on CRAN, click the green badge above.

### Overview

`spatstat.random` supports

- generating random spatial patterns of points according to many simple rules
(complete spatial randomness, binomial process, random grid,
systematic random, stratified random, 
simple sequential inhibition, cell process),

- randomised alteration of patterns (thinning,
random shift, jittering),

- generating simulated realisations of spatial point processes
(Poisson processes, Matern inhibition models, Matern cluster processes,
Neyman-Scott cluster processes, log-Gaussian Cox processes,
product shot noise cluster processes, Gibbs point processes)

- generating simulated realisations of Gibbs point processes
(Metropolis-Hastings birth-death-shift algorithm;
perfect simulation/ dominated coupling from the past;
alternating Gibbs sampler)

- generating random spatial patterns of line segments

- generating random tessellations

- generating random images (random noise, random mosaics).

Exceptions:

- generation of determinantal point processes is provided in `spatstat.model`

- generation of quasi-random patterns is provided in `spatstat.geom`


### Installing the package

This repository contains the _development version_ of
`spatstat.random`. The easiest way to install the development version
is to start R and type

```R
repo <- c('https://spatstat.r-universe.dev', 'https://cloud.r-project.org')
install.packages("spatstat.random", dependencies=TRUE, repos=repo)
```

To install the latest _public release_ of `spatstat.random`,
type

```R
install.packages("spatstat.random")
```

