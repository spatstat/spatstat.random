\name{runifdisc}
\alias{runifdisc}
\title{Generate N Uniform Random Points in a Disc}
\description{
  Generate a random point pattern
  containing \eqn{n} independent uniform random points
  in a circular disc.
}
\usage{
 runifdisc(n, radius=1, centre=c(0,0), ..., nsim=1, drop=TRUE)
}
\arguments{
  \item{n}{
    Number of points.
  }
  \item{radius}{Radius of the circle.}
  \item{centre}{Coordinates of the centre of the circle.}
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{disc}} controlling the
    accuracy of approximation to the circle.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"})
  if \code{nsim=1}, or a list of point patterns if \code{nsim > 1}.
}
\details{
  This function generates \code{n} independent random points,
  uniformly distributed in a circular disc.

  It is faster (for a circular window) than the general
  code used in \code{\link{runifpoint}}.

  To generate random points in an ellipse, first generate points in a
  circle using \code{runifdisc},
  then transform to an ellipse using \code{\link[spatstat.geom]{affine}},
  as shown in the examples.
  
  To generate random points in other windows, use
  \code{\link{runifpoint}}.
  To generate non-uniform random points, use \code{\link{rpoint}}.
}
\seealso{
\code{\link[spatstat.geom]{disc}},
\code{\link{runifpoint}},
\code{\link{rpoint}}
}
\examples{
 # 100 random points in the unit disc
 plot(runifdisc(100))
 # 42 random points in the ellipse with major axis 3 and minor axis 1
 X <- runifdisc(42)
 Y <- affine(X, mat=diag(c(3,1)))
 plot(Y)
}
\author{
  \adrian
  and \rolf.
}
\keyword{spatial}
\keyword{datagen}

