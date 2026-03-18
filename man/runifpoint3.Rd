\name{runifpoint3}
\alias{runifpoint3}
\title{
  Generate N Uniform Random Points in Three Dimensions
}
\description{
  Generate a random point pattern containing \code{n}
  independent, uniform random points in three dimensions.
}
\usage{
runifpoint3(n, domain = box3(), \dots, nsim=1, drop=TRUE, ex=NULL)
}
\arguments{
  \item{n}{
    Number of points to be generated.
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
    which will be used to determine the number of points
    \code{n} and the containing box \code{domain}.
  }
}
\value{
  If \code{nsim = 1} and \code{drop=TRUE}, a point pattern in
  three dimensions (an object of class \code{"pp3"}).
  If \code{nsim > 1}, a list of such point patterns.
}
\details{
  This function generates \code{n} independent random points,
  uniformly distributed in the three-dimensional box \code{domain}.
}
\seealso{
  \code{\link{rpoispp3}}, 
  \code{\link[spatstat.geom]{pp3}}, 
  \code{\link[spatstat.geom]{box3}}
}
\examples{
   X <- runifpoint3(50)

   runifpoint3(ex=osteo$pts[[1]])
}
\author{
  \adrian
  and \rolf.
}
\keyword{spatial}
\keyword{datagen}
\concept{Three-dimensional}
