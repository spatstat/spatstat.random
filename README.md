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
(Poisson processes, \Matern inhibition models, \Matern cluster processes,
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

### Detailed contents

#### Generating random patterns

- binomial random patterns (`runifpoint`, `rpoint`, `rmpoint`, `runifdisc`)
- completely random patterns (`rpoispp`, `rmpoispp`)
- systematic random patterns (`rstrat`, `rsyst`)

#### Point process models defined by simple constructions

- simple sequential inhibition (`rSSI`)
- \Matern inhibition models (`rMaternI`, `rMaternII`)
- cell process (`rcell`)

#### Randomly changing an existing point pattern

- random shift (`rshift`)
- random thinning (`rthin`)
- random (re)labelling (`rlabel`)
- block resampling (`quadratresample`)

#### Clustered point processes

- log-Gaussian Cox process (`rLGCP`)
- Neyman-Scott cluster processes
(`rThomas`, `rMatClust`, `rCauchy`, `rVarGamma`)
- general Neyman-Scott cluster process (`rNeymanScott`)
- general Poisson cluster process (`rPoissonCluster`)
- Gauss-Poisson process (`rGaussPoisson`)

#### Gibbs point processes

- perfect simulation algorithms for specific Gibbs models
(`rHardcore`, `rStrauss`, `rStraussHard`, `rDiggleGratton`, `rDGS`,
`rPenttinen`, 
- Metropolis-Hastings simulation algorithm for Gibbs models
(`rmh`)
- alternating Gibbs sampler for multitype Gibbs processes (`rags`,
`ragsMultiHard`)
- alternating Gibbs sampler for area-interaction process (`ragsAreaInter`)

#### random points along lines

- random points along specified line segments
(`runifpointOnLines`, `rpoisppOnLines`)

#### random pixel images and random sets

- random pixel noise (`rnoise`)
- random mosaic (`rMosaicField`, `rMosaicSet`)

#### random line segment patterns

- Poisson line process (`rpoisline`)

#### random tessellations

- tessellation using Poisson line process (`rpoislinetess`)

#### three-dimensional point patterns

- uniform random points in 3D (`runifpoint3`)
- Poisson point process in 3D (`rpoispp3`)

#### multi-dimensional point patterns

- uniform random points in space or space-time (`runifpointx`)
- Poisson point process in space or space-time (`rpoisppx`)

#### probability distributions

- theoretical distribution of nearest neighbour distance (`rknn`)
- mixed Poisson distribution (`dmixpois`)

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

