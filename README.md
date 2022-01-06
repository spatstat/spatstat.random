# spatstat.random

## Random Generation and Simulation for the spatstat family

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.random)](http://cran.r-project.org/web/packages/spatstat.random) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.random)](https://github.com/spatstat/spatstat.random)

The original `spatstat` package has been split into
several sub-packages.

This package `spatstat.random` is one of these packages.
It contains the functions for

- generating random spatial patterns of points according to many simple rules
(complete spatial randomness, Poisson, binomial, random grid,
systematic, cell),

- randomised alteration of patterns (thinning,
random shift, jittering),

- simulated realisations of random point processes
(simple sequential inhibition, Matern inhibition models, Matern cluster process,
Neyman-Scott cluster processes, log-Gaussian Cox processes,
product shot noise cluster processes)

- simulation of Gibbs point processes
(Metropolis-Hastings birth-death-shift algorithm;
perfect simulation/ dominated coupling from the past;
alternating Gibbs sampler)

- generating random spatial patterns of line segments

- generating random tessellations

- generating random images (random noise, random mosaics).

The reorganisation of `spatstat` into a family of packages is described
on the GitHub repository
[spatstat/spatstat](https://github.com/spatstat/spatstat).
