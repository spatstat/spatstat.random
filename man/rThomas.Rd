\name{rThomas}
\alias{rThomas}
\title{Simulate Thomas Process}
\description{
  Generate a random point pattern, a simulated realisation of the
  Thomas cluster process.
}
\usage{
  rThomas(kappa, scale, mu, win = square(1),
             nsim=1, drop=TRUE,
             \dots,
             algorithm=c("BKBC", "naive"),
             nonempty=TRUE, 
             poisthresh=1e-6,
             expand = 4*scale,
             saveparents=FALSE, saveLambda=FALSE,
             kappamax=NULL, mumax=NULL, sigma) 
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a function, or a pixel image.
  }
  \item{scale}{
    Cluster size.
    Standard deviation of random displacement (along each coordinate axis)
    of a point from its cluster centre. A single positive number.
  }
  \item{mu}{
    Mean number of points per cluster (a single positive number)
    or reference intensity for the cluster points (a function or
    a pixel image).
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link[spatstat.geom]{as.owin}}.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
  \item{\dots}{
    Passed to \code{\link[spatstat.random]{clusterfield}} to control the image
    resolution when \code{saveLambda=TRUE}.
  }
  \item{algorithm}{
    String (partially matched) specifying the simulation algorithm.
    See Details.
  }
  \item{nonempty}{
    Logical. If \code{TRUE} (the default), a more efficient algorithm is
    used, in which parents are generated conditionally on having at
    least one offspring point in the window.
    If \code{FALSE}, parents are generated
    even if they have no offspring in the window. The default
    is recommended unless you need to simulate all the parent points
    for some other purpose.
  }
  \item{poisthresh}{
    Numerical threshold below which the model will be treated
    as a Poisson process. See Details.
  }
  \item{expand}{
    Window expansion distance. A single number.
    The distance by which the original window will be expanded
    in order to generate parent points.
    Has a sensible default.
  }
  \item{saveparents}{
    Logical value indicating whether to save the locations of the
    parent points as an attribute.
  }
  \item{saveLambda}{
    Logical. If \code{TRUE} then the random intensity corresponding to
    the simulated parent points will also be calculated and saved,
    and returns as an attribute of the point pattern.
  }
  \item{kappamax}{
    Optional. Numerical value which is an upper bound for the
    values of \code{kappa}, when \code{kappa} is a pixel image or a
    function.
  }
  \item{mumax}{
    Optional. Numerical value which is an upper bound for the
    values of \code{mu}, when \code{mu} is a pixel image or a
    function.
  }
  \item{sigma}{
    Deprecated. Equivalent to \code{scale}.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"}) if \code{nsim=1},
  or a list of point patterns if \code{nsim > 1}.

  Additionally, some intermediate results of the simulation are returned
  as attributes of this point pattern (see \code{\link[spatstat.random]{rNeymanScott}}).
  Furthermore, the simulated intensity
  function is returned as an attribute \code{"Lambda"}, if
  \code{saveLambda=TRUE}.
}
\details{
  This algorithm generates a realisation of the (`modified')
  Thomas process, a special case of the Neyman-Scott process,
  inside the window \code{win}.

  In the simplest case, where \code{kappa} and \code{mu}
  are single numbers, the cluster process is formed by first
  generating a uniform Poisson point process of \dQuote{parent} points 
  with intensity \code{kappa}. Then each parent point is
  replaced by a random cluster of \dQuote{offspring} points,
  the number of points per cluster being Poisson (\code{mu})
  distributed, and their
  positions being isotropic Gaussian displacements from the
  cluster parent location. The resulting point pattern
  is a realisation of the classical
  \dQuote{stationary Thomas process} generated inside the window \code{win}.
  This point process has intensity \code{kappa * mu}.

  Note that, for correct simulation of the model,
  the parent points are not restricted to lie inside the 
  window \code{win};
  the parent process is effectively the uniform Poisson process
  on the infinite plane.

  The algorithm can also generate spatially inhomogeneous versions of
  the Thomas process:
  \itemize{
    \item The parent points can be spatially inhomogeneous.
    If the argument \code{kappa} is a \code{function(x,y)}
    or a pixel image (object of class \code{"im"}), then it is taken
    as specifying the intensity function of an inhomogeneous Poisson
    process that generates the parent points.
    \item The offspring points can be inhomogeneous. If the
    argument \code{mu} is a \code{function(x,y)}
    or a pixel image (object of class \code{"im"}), then it is
    interpreted as the reference density for offspring points,
    in the sense of Waagepetersen (2007).
    For a given parent point, the offspring constitute a Poisson process
    with intensity function equal to \code{mu * f},
    where \code{f} is the Gaussian probability density
    centred at the parent point. Equivalently we first generate,
    for each parent point, a Poisson (\code{mumax}) random number of
    offspring (where \eqn{M} is the maximum value of \code{mu})
    with independent Gaussian displacements from the parent
    location, and then randomly thin the offspring points, with
    retention probability \code{mu/M}.
    \item Both the parent points and the offspring points can be
    spatially inhomogeneous, as described above.
  }

  Note that if \code{kappa} is a pixel image, its domain must be larger
  than the window \code{win}. This is because an offspring point inside
  \code{win} could have its parent point lying outside \code{win}.
  In order to allow this, the simulation algorithm
  first expands the original window \code{win}
  by a distance \code{expand} and generates the Poisson process of
  parent points on this larger window. If \code{kappa} is a pixel image,
  its domain must contain this larger window.

  The intensity of the Thomas process is \code{kappa * mu}
  if either \code{kappa} or \code{mu} is a single number. In the general
  case the intensity is an integral involving \code{kappa}, \code{mu}
  and \code{f}.

  If the pair correlation function of the model is very close
  to that of a Poisson process, deviating by less than
  \code{poisthresh}, then the model is approximately a Poisson process,
  and will be simulated as a Poisson process with intensity
  \code{kappa * mu}, using \code{\link[spatstat.random]{rpoispp}}. 
  This avoids computations that would otherwise require huge amounts
  of memory.
}
\section{Simulation Algorithm}{
  Two simulation algorithms are implemented.
  \itemize{
    \item The \emph{naive} algorithm generates the cluster process 
    by directly following the description given above. First the window
    \code{win} is expanded by a distance equal to \code{expand}.
    Then the parent points are generated in the expanded window according to
    a Poisson process with intensity \code{kappa}. Then each parent
    point is replaced by a finite cluster of offspring points as
    described above.
    The naive algorithm is used if \code{algorithm="naive"} or if
    \code{nonempty=FALSE}.
    \item The \emph{BKBC} algorithm, proposed by Baddeley and Chang
    (2023), is a modification of the algorithm of Brix and Kendall (2002). 
    Parents are generated in the infinite plane, subject to the
    condition that they have at least one offspring point inside the
    window \code{win}.
    The BKBC algorithm is used when \code{algorithm="BKBC"} (the default)
    and \code{nonempty=TRUE} (the default).
  }
  The naive algorithm becomes very slow when \code{scale} is large,
  while the BKBC algorithm is uniformly fast (Baddeley and Chang, 2023).

  If \code{saveparents=TRUE}, then the simulated point pattern will
  have an attribute \code{"parents"} containing the coordinates of the
  parent points, and an attribute \code{"parentid"} mapping each
  offspring point to its parent.

  If \code{nonempty=TRUE} (the default), then parents are generated
  subject to the condition that they have at least one offspring point 
  in the window \code{win}. 
  \code{nonempty=FALSE}, then parents without offspring will be included;
  this option is not available in the \emph{BKBC} algorithm.
  
  Note that if \code{kappa} is a pixel image, its domain must be larger
  than the window \code{win}. This is because an offspring point inside
  \code{win} could have its parent point lying outside \code{win}.
  In order to allow this, the naive simulation algorithm
  first expands the original window \code{win}
  by a distance equal to \code{expand} and generates the Poisson process of
  parent points on this larger window. If \code{kappa} is a pixel image,
  its domain must contain this larger window.

  If the pair correlation function of the model is very close
  to that of a Poisson process, with maximum deviation less than
  \code{poisthresh}, then the model is approximately a Poisson process.
  This is detected by the naive algorithm which then
  simulates a Poisson process with intensity
  \code{kappa * mu}, using \code{\link[spatstat.random]{rpoispp}}. 
  This avoids computations that would otherwise require huge amounts
  of memory.
}
\section{Fitting cluster models to data}{
  The Thomas model with homogeneous parents
  (i.e. where \code{kappa} is a single number)
  where the offspring are either homogeneous or inhomogeneous (\code{mu}
  is a single number, a function or pixel image)
  can be fitted to point pattern data using \code{\link[spatstat.model]{kppm}},
  or fitted to the inhomogeneous \eqn{K} function
  using \code{\link[spatstat.model]{thomas.estK}}
  or \code{\link[spatstat.model]{thomas.estpcf}}.
  
  Currently \pkg{spatstat} does not support fitting the
  Thomas cluster process model
  with inhomogeneous parents.

  A Thomas cluster process model fitted by \code{\link[spatstat.model]{kppm}}
  can be simulated automatically using \code{\link[spatstat.model]{simulate.kppm}}
  (which invokes \code{rThomas} to perform the simulation).
}
\seealso{
\code{\link[spatstat.random]{rpoispp}},    
\code{\link[spatstat.random]{rMatClust}},  
\code{\link[spatstat.random]{rCauchy}},    
\code{\link[spatstat.random]{rVarGamma}},   
\code{\link[spatstat.random]{rNeymanScott}},  
\code{\link[spatstat.random]{rGaussPoisson}}.  

For fitting the model, see
\code{\link[spatstat.model]{kppm}},  
\code{\link[spatstat.model]{clusterfit}}.  
}

\references{
  \baddchangclustersim

  Brix, A. and Kendall, W.S. (2002)
  Simulation of cluster point processes without edge effects.
  \emph{Advances in Applied Probability} \bold{34}, 267--280.

  Diggle, P. J., Besag, J. and Gleaves, J. T. (1976)
  Statistical analysis of spatial point patterns by
  means of distance methods. \emph{Biometrics} \bold{32} 659--667.

  Thomas, M. (1949) A generalisation of Poisson's binomial limit for use
  in ecology. \emph{Biometrika} \bold{36}, 18--25.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\examples{
  #homogeneous
  X <- rThomas(10, 0.2, 5)
  #inhomogeneous
  Z <- as.im(function(x,y){ 5 * exp(2 * x - 1) }, owin())
  Y <- rThomas(10, 0.2, Z)
}
\author{
  \adrian, \rolf and \yamei.
}
\keyword{spatial}
\keyword{datagen}

