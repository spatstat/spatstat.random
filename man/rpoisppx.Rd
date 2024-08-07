\name{rpoisppx}
\alias{rpoisppx}
\title{
  Generate Poisson Point Pattern in Any Dimensions
}
\description{
  Generate a random multi-dimensional point pattern
  using the homogeneous Poisson process.
}
\usage{
rpoisppx(lambda, domain, nsim=1, drop=TRUE)
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process.
    A single positive number.
  }
  \item{domain}{
    Multi-dimensional box in which the process should be generated.
    An object of class \code{"boxx"}.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a single point pattern.
  }
}
\value{
  If \code{nsim = 1} and \code{drop=TRUE},
  a point pattern (an object of class \code{"ppx"}).
  If \code{nsim > 1} or \code{drop=FALSE}, a list of such point patterns.
}
\details{
  This function generates a realisation
  of the homogeneous Poisson process in multi dimensions,
  with intensity \code{lambda} (points per unit volume).
  
  The realisation is generated inside the multi-dimensional region
  \code{domain} which currently must be a rectangular box (object of
  class \code{"boxx"}).
}
\note{
  The intensity \code{lambda} is the expected number of points
  \emph{per unit volume}. 
}
\seealso{
  \code{\link{runifpointx}}, 
  \code{\link[spatstat.geom]{ppx}}, 
  \code{\link[spatstat.geom]{boxx}}
}
\examples{
   w <- boxx(x=c(0,1), y=c(0,1), z=c(0,1), t=c(0,3))
   X <- rpoisppx(10, w)
}
\author{
  \adrian
  and \rolf.
}
\keyword{spatial}
\keyword{datagen}
