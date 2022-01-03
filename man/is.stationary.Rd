\name{is.stationary}
\alias{is.stationary}
\alias{is.stationary.rmhmodel}
\alias{is.poisson}
\alias{is.poisson.rmhmodel}
\title{
  Recognise Stationary and Poisson Point Process Models
}
\description{
  Given a point process model (either a model that has been fitted to
  data, or a model specified by its parameters), 
  determine whether the model is a stationary point process,
  and whether it is a Poisson point process.
}
\usage{
is.stationary(x)
\method{is.stationary}{rmhmodel}(x)

is.poisson(x)
\method{is.poisson}{rmhmodel}(x)
}
\arguments{
  \item{x}{
    A fitted spatial point process model
    (object of class \code{"ppm"}, \code{"kppm"},
    \code{"lppm"}, \code{"dppm"} or \code{"slrm"})
    or a specification of a Gibbs point process model
    (object of class \code{"rmhmodel"})
    or a similar object.
  }
}
\details{
  The argument \code{x} represents a fitted spatial point process model
  or a similar object.

  \code{is.stationary(x)} returns \code{TRUE} if \code{x} represents
  a stationary point process, and \code{FALSE} if not.

  \code{is.poisson(x)} returns \code{TRUE} if \code{x} represents
  a Poisson point process, and \code{FALSE} if not.

  The functions \code{is.stationary} and \code{is.poisson} are generic,
  with methods for the classes \code{"ppm"} (Gibbs point process models),
  \code{"kppm"} (cluster or Cox point process models),
  \code{"slrm"} (spatial logistic regression models) and
  \code{"rmhmodel"} (model specifications for the
  Metropolis-Hastings algorithm).
  Additionally \code{is.stationary} has a method for
  classes \code{"detpointprocfamily"} and \code{"dppm"}
  (both determinantal point processes) and
  \code{is.poisson} has a method for 
  class \code{"interact"} (interaction structures for Gibbs models). 

  \code{is.poisson.kppm} will return \code{FALSE}, unless
  the model \code{x} is degenerate:
  either \code{x} has zero intensity so that its realisations are empty
  with probability 1, or it is a log-Gaussian Cox process
  where the log intensity has zero variance.

  \code{is.poisson.slrm} will always return \code{TRUE},
  by convention.
}
\value{
  A logical value.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{is.marked}} to determine whether a model is a marked
  point process. 

}
\examples{
  m <- rmhmodel(cif='strauss', par=list(beta=10, gamma=0.1, r=1))
  is.stationary(m)
  is.poisson(m)
  is.poisson(rmhmodel(cif='strauss', par=list(beta=10, gamma=1, r=1)))
}
\keyword{spatial}
\keyword{models}
