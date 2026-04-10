\name{rUnround}
\alias{rUnround}
\alias{rUnround.ppp}
\title{
  Random Un-Rounding of Spatial Location
}
\description{
  Given a point pattern in which the spatial coordinates have been
  discretised (rounded), randomly displace each point to reverse the
  effect of rounding.
}
\usage{
rUnround(X, \dots)

\method{rUnround}{ppp}(X, \dots, xstep=NULL, ystep=xstep, 
    nsim = 1, drop=TRUE, giveup = 1000)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} specifying
    the discretisation, if \code{xstep} and \code{ystep} are not given.
  }
  \item{xstep, ystep}{
    Numeric values giving the resolution of the discretised coordinates.
  }
  \item{giveup}{
    Maximum number of attempts to place a point inside the window.
  }
  \item{nsim}{
    Number of simulated realisations of the unrounded pattern.
  }
  \item{drop}{
    Logical value specifying, when \code{nsim=1}, whether to return a
    point pattern (\code{drop=TRUE}, the default) or a list containing
    one point pattern (\code{drop=FALSE}).
  }
}
\details{
  This function is similar to \code{\link[spatstat.geom]{rjitter}} except that it is
  designed to reverse the effect of discretisation of coordinates.
  
  The spatial coordinates of \code{X} are modified by
  adding independent, uniformly distributed random numbers
  to the \eqn{x} and \eqn{y} coordinates. If the points of \code{X}
  were all located at the centres of pixels in a pixel grid, then the points of
  \code{rUnround(X)} will be uniformly randomly distributed
  within the same pixels.
}
\value{
  Another point pattern (object of class \code{"ppp"})
  in the same window as \code{X} and with
  the same number of points as \code{X}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.geom]{rjitter}}, \code{\link[spatstat.geom]{discretise}}
}
\examples{
  X <- runifrect(30, Frame(letterR))[letterR]
  Y <- discretise(X, eps=0.05)
  Z <- rUnround(Y)
}
\keyword{spatial}
\keyword{datagen}
\keyword{manip}
