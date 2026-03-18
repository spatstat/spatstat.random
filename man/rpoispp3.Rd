\name{rpoispp3}
\alias{rpoispp3}
\title{
  Generate Poisson Point Pattern in Three Dimensions
}
\description{
  Generate a random three-dimensional point pattern
  using the homogeneous Poisson process.
}
\usage{
rpoispp3(lambda, domain = box3(), \dots, nsim=1, drop=TRUE, ex=NULL, lmax=NULL)
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process.
    A single positive number, or a \code{function(x,y,z)}.
  }
  \item{domain}{
    Three-dimensional box in which the process should be generated.
    An object of class \code{"box3"}.
  }
  \item{\dots}{Ignored.}
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
  \item{ex}{
    An example point pattern (object of class \code{"pp3"})
    which will be used to determine the intensity
    \code{lambda} and the containing box \code{domain}.
  }
  \item{lmax}{
    Maximum possible value of the intensity function \code{lambda},
    when \code{lambda} is a function.
  }
}
\value{
  If \code{nsim = 1} and \code{drop=TRUE}, a point pattern in
  three dimensions (an object of class \code{"pp3"}).
  If \code{nsim > 1}, a list of such point patterns.
}
\details{
  This function generates a realisation
  of the Poisson process in three dimensions,
  with intensity \code{lambda} (points per unit volume).

  If \code{lambda} is a single number, the homogeneous Poisson process
  with constant intensity \code{lambda} is generated.

  If \code{lambda} is a function with arguments \code{x,y,z}, the
  inhomogeneous Poisson process with intensity function
  \code{lambda(x,y,z)} is generated. In this case the argument
  \code{lmax} is required.
  
  The realisation is generated inside the three-dimensional region
  \code{domain} which currently must be a rectangular box (object of
  class \code{"box3"}).
}
\note{
  The intensity \code{lambda} is the expected number of points
  \emph{per unit volume}. 
}
\seealso{
  \code{\link{runifpoint3}}, 
  \code{\link[spatstat.geom]{pp3}}, 
  \code{\link[spatstat.geom]{box3}}
}
\examples{
   X <- rpoispp3(50)

   rpoispp3(ex=osteo$pts[[1]])
}
\author{
  \adrian
  and \rolf.
}
\keyword{spatial}
\keyword{datagen}
\concept{Three-dimensional}
