\name{rMosaicField}
\alias{rMosaicField}
\title{Mosaic Random Field}
\description{
  Generate a realisation of a random field
  which is piecewise constant on the tiles of a given tessellation.
}
\usage{
rMosaicField(X, 
    rgen = function(n) { sample(0:1, n, replace = TRUE)},
    ...,
    rgenargs=NULL)
}
\arguments{
  \item{X}{
    A tessellation (object of class \code{"tess"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}}
    determining the pixel resolution.
  }
  \item{rgen}{
    Function that generates random values for the
    tiles of the tessellation.
  }
  \item{rgenargs}{
    List containing extra arguments that should be passed
    to \code{rgen} (typically specifying parameters of the
    distribution of the values).
  }
}
\details{
  This function generates a realisation of a random field
  which is piecewise constant on the tiles of the given
  tessellation \code{X}. The values in each tile
  are independent and identically distributed.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\author{\adrian and \rolf}
\seealso{
  \code{\link{rpoislinetess}},
  \code{\link{rMosaicSet}}
}
\examples{
   if(interactive()) {
     lambda <- 3
     d <- 256
     n <- 30
   } else {
     lambda <- 1
     d <- 32
     n <- 5
   }
   X <- rpoislinetess(lambda)
   plot(rMosaicField(X, runif, dimyx=d))
   plot(rMosaicField(X, rnorm, rgenargs=list(mean=10, sd=2), dimyx=d))
   Y <- dirichlet(runifpoint(n))
   plot(rMosaicField(Y, rnorm, dimyx=d))
}
\keyword{spatial}
\keyword{datagen}
