\name{rVarGamma}
\alias{rVarGamma}
\title{Simulate Neyman-Scott Point Process with Variance Gamma cluster kernel}
\description{
  Generate a random point pattern, a simulated realisation of the
  Neyman-Scott process with Variance Gamma (Bessel) cluster kernel.
}
\usage{
 rVarGamma(kappa, scale, mu, 
              nu,
              win = square(1),
              nsim=1, drop=TRUE,
              \dots,
              n.cond = NULL, w.cond = NULL,
              algorithm=c("BKBC", "naive"),
              nonempty=TRUE, 
              thresh = 0.001,
              poisthresh=1e-6,
              expand = NULL,
              saveparents=FALSE, saveLambda=FALSE,
              kappamax=NULL, mumax=NULL, LambdaOnly=FALSE) 
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a function, or a pixel image.
  }
  \item{scale}{
    Scale parameter for cluster kernel. Determines the size of clusters.
    A single positive number, in the same units as the spatial coordinates.
  }
  \item{mu}{
    Mean number of points per cluster (a single positive number)
    or reference intensity for the cluster points (a function or
    a pixel image).
  }
  \item{nu}{
    Shape parameter for the cluster kernel. A number greater than -1.
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
    resolution when \code{saveLambda=TRUE}, and to 
    \code{\link[spatstat.random]{clusterradius}} when
    \code{expand} is missing or \code{NULL}.
  }
  \item{n.cond}{
    Optional. Integer specifying a fixed number of points.
    See the section on \emph{Conditional Simulation}.
  }
  \item{w.cond}{
    Optional. Conditioning region. A window (object of class \code{"owin"})
    specifying the region which must contain exactly \code{n.cond} points.
    See the section on \emph{Conditional Simulation}.
  }
  \item{algorithm}{
    String (partially matched) specifying the simulation algorithm.
    See Details.
  }
  \item{nonempty}{
    Logical. If \code{TRUE} (the default), a more efficient algorithm is
    used, in which parents are generated conditionally on having at
    least one offspring point. If \code{FALSE}, parents are generated
    even if they have no offspring. Both choices are valid; the default
    is recommended unless you need to simulate all the parent points
    for some other purpose.
  }
  \item{thresh}{
    Threshold relative to the cluster kernel value at the origin (parent
    location) determining when the cluster kernel will be treated as
    zero for simulation purposes. Will be overridden by argument
    \code{expand} if that is given.
  }
  \item{poisthresh}{
    Numerical threshold below which the model will be treated
    as a Poisson process. See Details.
  }
  \item{expand}{
    Window expansion distance. A single number.
    The distance by which the original window will be expanded
    in order to generate parent points.
    Has a sensible default, determined by calling
    \code{\link[spatstat.random]{clusterradius}}
    with the numeric threshold value given
    in \code{thresh}.
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
  \item{LambdaOnly}{
    Logical value specifying whether to return only the random intensity,
    rather than the point pattern.
  }
}
\value{
  A point pattern (object of class \code{"ppp"})
  or a list of point patterns.

  Additionally, some intermediate results of the simulation are returned
  as attributes of this point pattern (see
  \code{\link[spatstat.random]{rNeymanScott}}). 
  Furthermore, the simulated intensity
  function is returned as an attribute \code{"Lambda"}, if
  \code{saveLambda=TRUE}.

  If \code{LambdaOnly=TRUE} the result is a pixel image
  (object of class \code{"im"}) or a list of pixel images.
}
\details{
  This algorithm generates a realisation of the Neyman-Scott process
  with Variance Gamma (Bessel) cluster kernel, inside the window \code{win}.

  The process is constructed by first
  generating a Poisson point process of ``parent'' points 
  with intensity \code{kappa}. Then each parent point is
  replaced by a random cluster of points, the number of points in each
  cluster being random with a Poisson (\code{mu}) distribution,
  and the points being placed independently and uniformly
  according to a Variance Gamma kernel.

  Note that, for correct simulation of the model,
  the parent points are not restricted to lie inside the 
  window \code{win};
  the parent process is effectively the uniform Poisson process
  on the infinite plane.

  The shape of the kernel is determined by the dimensionless
  index \code{nu}. This is the parameter \eqn{nu^\prime}{nu'} (nu-prime)
  \eqn{\nu^\prime = \alpha/2-1}{nu' = alpha/2 - 1} appearing in
  equation (12) on page 126 of Jalilian et al (2013). 

  The scale of the kernel is determined by the argument \code{scale},
  which is the parameter
  \eqn{\eta}{eta} appearing in equations (12) and (13) of
  Jalilian et al (2013). 
  It is expressed in units of length (the same as the unit of length for 
  the window \code{win}).
  
  The algorithm can also generate spatially inhomogeneous versions of
  the cluster process:
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
    in the sense of Waagepetersen (2006).
  }
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
\section{Conditional Simulation}{
  If \code{n.cond} is specified, it should be a single integer.
  Simulation will be conditional on the event 
  that the pattern contains exactly \code{n.cond} points
  (or contains exactly \code{n.cond} points inside
  the region \code{w.cond} if it is given).

  Conditional simulation uses the rejection algorithm described
  in Section 6.2 of Moller, Syversveen and Waagepetersen (1998).
  There is a maximum number of proposals which will be attempted.
  Consequently the return value may contain fewer
  than \code{nsim} point patterns.

  The current algorithm for conditional simulation ignores the
  argument \code{saveparents} and does not save the parent points.
}
\section{Fitting cluster models to data}{
  The Variance-Gamma cluster model with homogeneous parents
  (i.e. where \code{kappa} is a single number)
  where the offspring are either homogeneous or inhomogeneous (\code{mu}
  is a single number, a function or pixel image)
  can be fitted to point pattern data using \code{\link[spatstat.model]{kppm}},
  or fitted to the inhomogeneous \eqn{K} function
  using \code{\link[spatstat.model]{vargamma.estK}}
  or \code{\link[spatstat.model]{vargamma.estpcf}}.
  
  Currently \pkg{spatstat} does not support fitting the
  Variance-Gamma cluster process model
  with inhomogeneous parents.

  A Variance-Gamma cluster process model fitted by \code{\link[spatstat.model]{kppm}}
  can be simulated automatically using \code{\link[spatstat.model]{simulate.kppm}}
  (which invokes \code{rVarGamma} to perform the simulation).
}
\section{Warning}{
  The argument \code{nu} is the parameter \eqn{\nu^\prime}{nu'} (nu-prime)
  \eqn{\nu^\prime = \alpha/2-1}{nu' = alpha/2 - 1} appearing in
  equation (12) on page 126 of Jalilian et al (2013).
  This is different from the parameter called \eqn{\nu}{nu} appearing in
  equation(14) on page 127 of  Jalilian et al (2013), defined
  by \eqn{\nu = \alpha-1}{nu = alpha - 1}. This has been a frequent
  source of confusion.
}
\seealso{
\code{\link[spatstat.random]{rpoispp}},    
\code{\link[spatstat.random]{rMatClust}},  
\code{\link[spatstat.random]{rThomas}},    
\code{\link[spatstat.random]{rCauchy}},   
\code{\link[spatstat.random]{rNeymanScott}},  
\code{\link[spatstat.random]{rGaussPoisson}}.  

For fitting the model, see
\code{\link[spatstat.model]{kppm}},  
\code{\link[spatstat.model]{clusterfit}},  
\code{\link[spatstat.model]{vargamma.estK}}, 
\code{\link[spatstat.model]{vargamma.estpcf}}. 
}
\examples{
 # homogeneous
 X <- rVarGamma(kappa=5, scale=2, mu=5, nu=-1/4)
 # inhomogeneous
 ff <- function(x,y){ exp(2 - 3 * abs(x)) }
 fmax <- exp(2)
 Z <- as.im(ff, W= owin())
 Y <- rVarGamma(kappa=5, scale=2, mu=Z, nu=0)
 YY <- rVarGamma(kappa=ff, scale=2, mu=3, nu=0, kappamax=fmax)
}
\references{
  \baddchangclustersim

  Brix, A. and Kendall, W.S. (2002)
  Simulation of cluster point processes without edge effects.
  \emph{Advances in Applied Probability} \bold{34}, 267--280.

  Jalilian, A., Guan, Y. and Waagepetersen, R. (2013)
  Decomposition of variance for spatial Cox processes.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 119-137.
  
  \Moller, J., Syversveen, A. and Waagepetersen, R. (1998)
  Log Gaussian Cox Processes.
  \emph{Scandinavian Journal of Statistics} \bold{25}, 451--482.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\author{
  Original algorithm by Abdollah Jalilian and Rasmus Waagepetersen.
  Adapted for \pkg{spatstat} by \adrian.
  Brix-Kendall-Baddeley-Chang algorithm implemented by
  \adrian and \yamei.
}
\keyword{spatial}
\keyword{datagen}

