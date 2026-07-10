\name{rpoisDirichletTess}
\alias{rpoisDirichletTess}
\title{Poisson Dirichlet Tessellation}
\description{
  Divide a random tessellation of a given spatial region,
  using the Dirichlet tessellation defined by
  a Poisson point process on the infinite plane. 
}
\usage{
rpoisDirichletTess(lambda, win = owin())
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson point process.
    Either a single positive number, or a function which can be
    evaluated at any spatial location.
  }
  \item{win}{
    Window to be divided into tiles.
    An object of class \code{"owin"}
    or something acceptable to \code{\link[spatstat.geom]{as.owin}}.
  }
}
\details{
  This algorithm computes the Dirichlet-Voronoi tessellation
  of the window \code{win} induced by a Poisson point process
  on the infinite two-dimensional plane, avoiding edge effects.
  
  The argument \code{lambda} is the intensity of the Poisson process.
  It should be either a single positive number, or a
  \code{function(x,y)}
  which can be evaluated at any spatial location \eqn{(x,y)}.

  The algorithm first generates a realisation of the Poisson process
  inside \code{win}, and computes the Dirichlet-Voronoi tessellation
  of this point pattern. Then in order to avoid edge effects,
  the algorithm calculates the
  maximum possible distance \eqn{d} outside the window
  such that a random point outside \code{w} could affect the
  outcome inside \code{w}. The window is expanded by this distance
  \eqn{d}; the realisation of the Poisson process in the expanded region
  is generated; and the Dirichlet-Voronoi tessellation of the augmented
  point pattern is computed, then restricted to the original window.

  If \code{lambda} is a single number, the result is a realisation
  of a stationary random tessellation (restricted to \code{win}).
}
\value{
  A tessellation (object of class \code{"tess"}).

  Also has an attribute \code{"X"} containing the realisation of the
  Poisson point process (restricted to those points which affected the
  result).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{rpoislinetess}}
}
\examples{
 A <- rpoisDirichletTess(10)
 X <- attr(A, "X")
 plot(Frame(X), type="n", main="rpoisDirichletTess(10)")
 plot(A, add=TRUE)
 plot(X, add=TRUE, pch=16)
}
\keyword{spatial}
\keyword{datagen}
\concept{Tessellation}
