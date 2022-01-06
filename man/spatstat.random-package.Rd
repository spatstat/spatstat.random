\name{spatstat.random-package}
\alias{spatstat.random-package} 
\alias{spatstat.random} 
\docType{package}
\title{The spatstat.random Package}
\description{
  The \pkg{spatstat.random} package
  belongs to the \pkg{spatstat} family of packages.
  It contains the functionality for generating random spatial patterns
  and simulation of random point processes.
}
\details{
  \pkg{spatstat} is a family of \R packages
  for the statistical analysis of spatial data.
  Its main focus is the analysis of
  spatial patterns of points in two-dimensional space.

  This sub-package \pkg{spatstat.random} contains the
  functions that perform random generation of spatial patterns
  and simulation of random point processes:
  \itemize{
    \item
    generating random spatial patterns of points
    according to many simple rules
    (complete spatial randomness, Poisson, binomial, random grid,
    systematic, cell);
    \item
    randomised alteration of spatial patterns (thinning, random shifting,
    jittering, random labelling);
    \item
    simulated realisations of random point processes
    (simple sequential inhibition,
    \Matern inhibition models,
    \Matern cluster process,
    Neyman-Scott cluster processes,
    log-Gaussian Cox processes,
    product shot noise cluster processes);
    \item
    simulation of Gibbs point processes
    (Metropolis-Hastings birth-death-shift algorithm,
    perfect simulation/ dominated coupling from the past,
    alternating Gibbs sampler);
    \item
    generating random spatial patterns of line segments;
    \item
    generating random tessellations;
    \item
    generating random images (random noise, random mosaics).
  }
  (Functions for linear networks
  are in the separate sub-package \pkg{spatstat.linnet}.)
}
\section{Structure of the spatstat family}{  
  The \pkg{spatstat} family of packages
  is designed to support a complete statistical analysis
  of spatial data. It supports
  \itemize{
    \item creation, manipulation and plotting of point patterns;
    \item exploratory data analysis;
    \item spatial random sampling;
    \item simulation of point process models;
    \item parametric model-fitting;
    \item non-parametric smoothing and regression;
    \item formal inference (hypothesis tests, confidence intervals);
    \item model diagnostics.
  }
  The orginal \pkg{spatstat} package grew to be very large.
  It has now been divided into several \bold{sub-packages}:
  \itemize{
    \item \pkg{spatstat.utils} containing basic utilities
    \item \pkg{spatstat.sparse} containing linear algebra utilities
    \item \pkg{spatstat.data} containing datasets
    \item \pkg{spatstat.geom} containing geometrical objects
    and geometrical operations
    \item \pkg{spatstat.random} containing functionality for
    simulation and random generation
    \item \pkg{spatstat.core} containing the main functionality
    for statistical analysis and modelling of spatial data
    \item \pkg{spatstat.linnet} containing functions for
    spatial data on a linear network
    \item \pkg{spatstat}, which simply loads the other sub-packages
    listed above, and provides documentation.
  }
  When you install \pkg{spatstat}, these sub-packages are also
  installed. Then if you load the \pkg{spatstat} package by typing
  \code{library(spatstat)}, the other sub-packages listed above will
  automatically be loaded or imported.

  For an overview of all the functions available in
  the sub-packages of \pkg{spatstat},
  see the help file for \code{"spatstat-package"}
  in the \pkg{spatstat} package.

  Additionally there are several \bold{extension packages:}
  \itemize{
    \item \pkg{spatstat.gui} for interactive graphics
    \item \pkg{spatstat.local} for local likelihood
    (including geographically weighted regression)
    \item \pkg{spatstat.Knet} for additional, computationally efficient code
    for linear networks
    \item \pkg{spatstat.sphere} (under development) for spatial data
    on a sphere, including spatial data on the earth's surface
  }
  The extension packages must be installed separately
  and loaded explicitly if needed. They also have separate documentation.
}
\section{Functionality in \pkg{spatstat.random}}{
  Following is a list of the functionality provided in the
  \pkg{spatstat.random} package only.
  
  \bold{To simulate a random point pattern:}
  
  \tabular{ll}{
	\code{\link[spatstat.random]{runifpoint}} \tab
	generate \eqn{n} independent uniform random points \cr
	\code{\link[spatstat.random]{rpoint}} \tab
	generate \eqn{n} independent random points \cr
	\code{\link[spatstat.random]{rmpoint}} \tab
	generate \eqn{n} independent multitype random points \cr
	\code{\link[spatstat.random]{rpoispp}} \tab
	simulate the (in)homogeneous Poisson point process \cr
	\code{\link[spatstat.random]{rmpoispp}} \tab
	simulate the (in)homogeneous multitype Poisson point process \cr
	\code{\link[spatstat.random]{runifdisc}} \tab
	generate \eqn{n} independent uniform random points in disc\cr
	\code{\link[spatstat.random]{rstrat}} \tab
	stratified random sample of points \cr
	\code{\link[spatstat.random]{rMaternI}}  \tab
	simulate the \Matern Model I inhibition process\cr
	\code{\link[spatstat.random]{rMaternII}} \tab
	simulate the \Matern Model II inhibition process\cr
	\code{\link[spatstat.random]{rSSI}} \tab
	simulate Simple Sequential Inhibition process\cr
	\code{\link[spatstat.random]{rStrauss}} \tab
	simulate Strauss process (perfect simulation)\cr
	\code{\link[spatstat.random]{rHardcore}} \tab
	simulate Hard Core process (perfect simulation)\cr
	\code{\link[spatstat.random]{rStraussHard}} \tab
	simulate Strauss-hard core process (perfect simulation)\cr
	\code{\link[spatstat.random]{rDiggleGratton}} \tab
	simulate Diggle-Gratton process (perfect simulation)\cr
	\code{\link[spatstat.random]{rDGS}} \tab
	simulate Diggle-Gates-Stibbard process (perfect simulation)\cr
	\code{\link[spatstat.random]{rPenttinen}} \tab
	simulate Penttinen process (perfect simulation)\cr
	\code{\link[spatstat.random]{rNeymanScott}} \tab
	simulate a general Neyman-Scott process\cr
	\code{\link[spatstat.random]{rPoissonCluster}} \tab
	simulate a general Poisson cluster process\cr
	\code{\link[spatstat.random]{rMatClust}} \tab
	simulate the \Matern Cluster process\cr
	\code{\link[spatstat.random]{rThomas}} \tab
	simulate the Thomas process  \cr
	\code{\link[spatstat.random]{rGaussPoisson}}  \tab
	simulate the Gauss-Poisson cluster process\cr
	\code{\link[spatstat.random]{rCauchy}} \tab
	simulate Neyman-Scott Cauchy cluster process \cr
	\code{\link[spatstat.random]{rVarGamma}} \tab
	simulate Neyman-Scott Variance Gamma cluster process \cr
	\code{\link[spatstat.random]{rthin}} \tab  random thinning  \cr
	\code{\link[spatstat.random]{rcell}} \tab
	simulate the Baddeley-Silverman cell process  \cr
	\code{\link[spatstat.random]{rmh}} \tab
	simulate Gibbs point process using Metropolis-Hastings \cr
	\code{\link[spatstat.random]{runifpointOnLines}} \tab
	generate \eqn{n} random points along specified line segments \cr
	\code{\link[spatstat.random]{rpoisppOnLines}} \tab
	generate Poisson random points along specified line segments 
      }

      \bold{To randomly change an existing point pattern:}
      
      \tabular{ll}{
	\code{\link[spatstat.random]{rshift}} \tab random shifting of points \cr
	\code{\link[spatstat.random]{rthin}} \tab  random thinning \cr
	\code{\link[spatstat.random]{rlabel}} \tab random (re)labelling of a multitype
	point pattern \cr
        \code{\link[spatstat.random]{quadratresample}} \tab block resampling 
      }

      \bold{Random pixel images:}

      An object of class \code{"im"} represents a pixel image.

      \tabular{ll}{
	\code{\link[spatstat.random]{rnoise}} \tab random pixel noise
     }

    \bold{Line segment patterns}

    An object of class \code{"psp"} represents a pattern of straight line
    segments.

    \tabular{ll}{
      \code{\link[spatstat.random]{rpoisline}} \tab generate a realisation of the
      Poisson line process inside a window
    }

    \bold{Tessellations}

    An object of class \code{"tess"} represents a tessellation.

    \tabular{ll}{
      \code{\link[spatstat.random]{rpoislinetess}} \tab generate tessellation using Poisson line
      process 
    }

    \bold{Three-dimensional point patterns}

    An object of class \code{"pp3"} represents a three-dimensional
    point pattern in a rectangular box. The box is represented by
    an object of class \code{"box3"}.

    \tabular{ll}{
      \code{\link[spatstat.random]{runifpoint3}} \tab generate uniform random points in 3-D \cr
      \code{\link[spatstat.random]{rpoispp3}} \tab generate Poisson random points in 3-D \cr
    }

    \bold{Multi-dimensional space-time point patterns}

    An object of class \code{"ppx"} represents a 
    point pattern in multi-dimensional space and/or time.

    \tabular{ll}{
      \code{\link[spatstat.random]{runifpointx}} \tab generate uniform random points \cr
      \code{\link[spatstat.random]{rpoisppx}} \tab generate Poisson random points 
    }
    

    \bold{Probability Distributions}

    \tabular{ll}{
      \code{\link[spatstat.random]{rknn}} \tab theoretical distribution of
      nearest neighbour distance \cr
      \code{\link[spatstat.random]{dmixpois}} \tab mixed Poisson distribution \cr
 }

  
  \bold{Simulation}

  There are many ways to generate a random point pattern,
  line segment pattern, pixel image or tessellation
  in \pkg{spatstat}. 

  \bold{Random point patterns:}

  \tabular{ll}{
    \code{\link[spatstat.random]{runifpoint}} \tab
    generate \eqn{n} independent uniform random points \cr
    \code{\link[spatstat.random]{rpoint}} \tab
    generate \eqn{n} independent random points \cr
    \code{\link[spatstat.random]{rmpoint}} \tab
    generate \eqn{n} independent multitype random points \cr
    \code{\link[spatstat.random]{rpoispp}} \tab
    simulate the (in)homogeneous Poisson point process \cr
    \code{\link[spatstat.random]{rmpoispp}} \tab
    simulate the (in)homogeneous multitype Poisson point process \cr
    \code{\link[spatstat.random]{runifdisc}} \tab
    generate \eqn{n} independent uniform random points in disc\cr
    \code{\link[spatstat.random]{rstrat}} \tab
    stratified random sample of points \cr
    \code{\link[spatstat.random]{rMaternI}}  \tab
    simulate the \Matern Model I inhibition process\cr
    \code{\link[spatstat.random]{rMaternII}} \tab
    simulate the \Matern Model II inhibition process\cr
    \code{\link[spatstat.random]{rSSI}} \tab
    simulate Simple Sequential Inhibition process\cr
    \code{\link[spatstat.random]{rHardcore}} \tab
    simulate hard core process (perfect simulation)\cr
    \code{\link[spatstat.random]{rStrauss}} \tab
    simulate Strauss process (perfect simulation)\cr
    \code{\link[spatstat.random]{rStraussHard}} \tab
    simulate Strauss-hard core process (perfect simulation)\cr
    \code{\link[spatstat.random]{rDiggleGratton}} \tab
    simulate Diggle-Gratton process (perfect simulation)\cr
    \code{\link[spatstat.random]{rDGS}} \tab
    simulate Diggle-Gates-Stibbard process (perfect simulation)\cr
    \code{\link[spatstat.random]{rPenttinen}} \tab
    simulate Penttinen process (perfect simulation)\cr
    \code{\link[spatstat.random]{rNeymanScott}} \tab
    simulate a general Neyman-Scott process\cr
    \code{\link[spatstat.random]{rMatClust}} \tab
    simulate the \Matern Cluster process\cr
    \code{\link[spatstat.random]{rThomas}} \tab
    simulate the Thomas process  \cr
    \code{\link[spatstat.random]{rLGCP}} \tab
    simulate the log-Gaussian Cox process  \cr
    \code{\link[spatstat.random]{rGaussPoisson}}  \tab
    simulate the Gauss-Poisson cluster process\cr
    \code{\link[spatstat.random]{rCauchy}} \tab
    simulate Neyman-Scott process with Cauchy clusters \cr
    \code{\link[spatstat.random]{rVarGamma}} \tab
    simulate Neyman-Scott process with Variance Gamma clusters \cr
    \code{\link[spatstat.random]{rcell}} \tab
    simulate the Baddeley-Silverman cell process  \cr
    \code{\link[spatstat.random]{runifpointOnLines}} \tab
    generate \eqn{n} random points along specified line segments \cr
    \code{\link[spatstat.random]{rpoisppOnLines}} \tab
    generate Poisson random points along specified line segments 
  }
      
  \bold{Resampling a point pattern:}

  \tabular{ll}{
    \code{\link[spatstat.random]{quadratresample}} \tab block resampling \cr
    \code{\link[spatstat.random]{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link[spatstat.random]{rthin}} \tab  random thinning 
  }
  
   \bold{Other random patterns:}

   \tabular{ll}{
     \code{\link[spatstat.random]{rpoisline}} \tab
     simulate the Poisson line process within a window \cr
     \code{\link[spatstat.random]{rpoislinetess}} \tab
     generate random tessellation using Poisson line process \cr
     \code{\link[spatstat.random]{rMosaicSet}} \tab
     generate random set by selecting some tiles of a tessellation \cr
     \code{\link[spatstat.random]{rMosaicField}} \tab
     generate random pixel image by assigning random values
     in each tile of a tessellation
   }

  \bold{Resampling and randomisation procedures}

  You can build your own tests based on randomisation
  and resampling using the following capabilities:
  
  \tabular{ll}{
    \code{\link[spatstat.random]{quadratresample}} \tab block resampling \cr
    \code{\link[spatstat.random]{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link[spatstat.random]{rthin}} \tab  random thinning  
  }
}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  \spatstatAuthors.
}
\section{Acknowledgements}{
  Kasper Klitgaard Berthelsen,
  Tilman Davies,
  Ute Hahn,
  Abdollah Jalilian,
  and
  Rasmus Waagepetersen
  made substantial contributions of code.

  For comments, corrections, bug alerts and suggestions, we thank
  Monsuru Adepeju,
  Corey Anderson,
  Ang Qi Wei,
  Ryan Arellano,
  Jens \ifelse{latex}{\out{{\AA}str{\" o}m}}{Astrom},
  Robert Aue,
  Marcel Austenfeld,
  Sandro Azaele,
  Malissa Baddeley,
  Guy Bayegnak,
  Colin Beale,
  Melanie Bell,
  Thomas Bendtsen,
  Ricardo Bernhardt,
  Andrew Bevan,
  Brad Biggerstaff,
  Anders Bilgrau,
  Leanne Bischof,
  Christophe Biscio,
  Roger Bivand,
  Jose M. Blanco Moreno,
  Florent Bonneu,
  Jordan Brown,
  Ian Buller,
  Julian Burgos,
  Simon Byers,
  Ya-Mei Chang,
  Jianbao Chen,
  Igor Chernayavsky,
  Y.C. Chin,
  Bjarke Christensen,
  \ifelse{latex}{\out{Luc\'{\i{}}a}}{Lucia} Cobo Sanchez,
  Jean-Francois Coeurjolly,
  Kim Colyvas,
  Hadrien Commenges,
  Rochelle Constantine,
  Robin Corria Ainslie,
  Richard Cotton,
  Marcelino de la Cruz,
  Peter Dalgaard,
  Mario D'Antuono,
  Sourav Das,
  Peter Diggle,
  Patrick Donnelly,
  Ian Dryden,
  Stephen Eglen,
  Ahmed El-Gabbas,
  Belarmain Fandohan,
  Olivier Flores,
  David Ford,
  Peter Forbes,
  Shane Frank,
  Janet Franklin,
  Funwi-Gabga Neba,
  Oscar Garcia,
  Agnes Gault,
  Jonas Geldmann,
  Marc Genton,
  Shaaban Ghalandarayeshi,
  Jason Goldstick,
  Pavel Grabarnik,
  C. Graf,
  Ute Hahn,
  Andrew Hardegen,
  Martin \Bogsted Hansen,
  Martin Hazelton,
  Juha Heikkinen,
  Mandy Hering,
  Markus Herrmann,
  Maximilian Hesselbarth,
  Paul Hewson,
  Hamidreza Heydarian,
  Kurt Hornik,
  Philipp Hunziker,
  Jack Hywood,
  Ross Ihaka,
  \ifelse{latex}{\out{\u{C}enk I\c{c}\"{o}s}}{Cenk Icos},
  Aruna Jammalamadaka,
  Robert John-Chandran,
  Devin Johnson,
  Mahdieh Khanmohammadi,
  Bob Klaver,
  Lily Kozmian-Ledward,
  Peter Kovesi,
  Mike Kuhn,
  Jeff Laake,
  Robert Lamb,
  \ifelse{latex}{\out{Fr\'{e}d\'{e}ric}}{Frederic} Lavancier,
  Tom Lawrence,
  Tomas Lazauskas,
  Jonathan Lee,
  George Leser,
  Angela Li,
  Li Haitao,
  George Limitsios,
  Andrew Lister,
  Nestor Luambua,
  Ben Madin,
  Martin Maechler,
  Kiran Marchikanti,
  Jeff Marcus,
  Robert Mark,
  Peter McCullagh,
  Monia Mahling,
  Jorge Mateu Mahiques,
  Ulf Mehlig,
  Frederico Mestre,
  Sebastian Wastl Meyer,
  Mi Xiangcheng,
  Lore De Middeleer,
  Robin Milne,
  Enrique Miranda,
  Jesper \Moller,
  Annie \ifelse{latex}{\out{Molli{\'e}}}{Mollie},
  Ines Moncada,
  Mehdi Moradi,
  Virginia Morera Pujol,
  Erika Mudrak,
  Gopalan Nair,
  Nader Najari,
  Nicoletta Nava,
  Linda Stougaard Nielsen,
  Felipe Nunes,
  Jens Randel Nyengaard,
  Jens \Oehlschlaegel,
  Thierry Onkelinx,
  Sean O'Riordan,
  Evgeni Parilov,
  Jeff Picka,
  Nicolas Picard,
  Tim Pollington,
  Mike Porter,
  Sergiy Protsiv,
  Adrian Raftery,
  Ben Ramage,
  Pablo Ramon,
  Xavier Raynaud,
  Nicholas Read,
  Matt Reiter,
  Ian Renner,
  Tom Richardson,
  Brian Ripley,
  Ted Rosenbaum,
  Barry Rowlingson,
  Jason Rudokas,
  Tyler Rudolph,
  John Rudge,
  Christopher Ryan,
  Farzaneh Safavimanesh,
  Aila \Sarkka,
  Cody Schank,
  Katja Schladitz,
  Sebastian Schutte,
  Bryan Scott,
  Olivia Semboli,
  \ifelse{latex}{\out{Fran\c{c}ois S\'{e}m\'{e}curbe}}{Francois Semecurbe},
  Vadim Shcherbakov,
  Shen Guochun,
  Shi Peijian,
  Harold-Jeffrey Ship,
  Tammy L Silva,
  Ida-Maria Sintorn,
  Yong Song, 
  Malte Spiess,
  Mark Stevenson,
  Kaspar Stucki,
  Jan Sulavik,
  Michael Sumner,
  P. Surovy,
  Ben Taylor,
  Thordis Linda Thorarinsdottir,
  Leigh Torres,
  Berwin Turlach,
  Torben Tvedebrink,
  Kevin Ummer,
  Medha Uppala,
  Andrew van Burgel,
  Tobias Verbeke,
  Mikko Vihtakari,
  Alexendre Villers,
  Fabrice Vinatier,
  Maximilian Vogtland,
  Sasha Voss,
  Sven Wagner,
  Hao Wang,
  H. Wendrock,
  Jan Wild,
  Carl G. Witthoft,
  Selene Wong,
  Maxime Woringer,
  Luke Yates,
  Mike Zamboni
  and
  Achim Zeileis.
}
\keyword{spatial}
\keyword{package}

