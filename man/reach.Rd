\name{reach}
\alias{reach}
\alias{reach.rmhmodel}
\title{Interaction Distance of a Point Process Model}
\description{
  Computes the interaction distance of a point process model.
}
\usage{
  reach(x, \dots)

  \method{reach}{rmhmodel}(x, \dots)

}
\arguments{
  \item{x}{Either a fitted point process model (object of class
    \code{"ppm"}), an interpoint interaction (object of class
    \code{"interact"}), a fitted interpoint interaction (object of
    class \code{"fii"}) or a point process model for simulation
    (object of class \code{"rmhmodel"}).
  }
  \item{\dots}{
    Other arguments are ignored.
  }
}
\value{
  The interaction distance, or \code{NA} if this cannot be
  computed from the information given.
}
\details{
  The function \code{reach} computes the 
  `interaction distance' or `interaction range' of a point process
  model.

  The definition of the interaction distance depends on the
  type of point process model. This help page explains the
  interaction distance for a Gibbs point process. For other kinds of
  models, see \code{\link[spatstat.model]{reach.kppm}} and
  \code{\link[spatstat.model]{reach.dppm}}.

  For a Gibbs point process model, the interaction distance
  is the shortest distance \eqn{D} such that any two points in the
  process which are separated by a distance greater than \eqn{D} do not
  interact with each other.

  For example, the interaction range of a Strauss process
  (see \code{\link[spatstat.model]{Strauss}} or \code{\link[spatstat.random]{rStrauss}})
  with parameters \eqn{\beta,\gamma,r}{beta,gamma,r} is equal to
  \eqn{r}, unless \eqn{\gamma=1}{gamma=1} in which case the model is
  Poisson and the interaction
  range is \eqn{0}.
  The interaction range of a Poisson process is zero.
  The interaction range of the Ord threshold process
  (see \code{\link[spatstat.model]{OrdThresh}}) is infinite, since two points \emph{may}
  interact at any distance apart.

  The function \code{reach} is generic, with methods
  for the case where \code{x} is 
  \itemize{
    \item
    a fitted point process model
    (object of class \code{"ppm"}, usually obtained from the model-fitting
    function \code{\link[spatstat.model]{ppm}});
    \item
    an interpoint interaction structure (object of class
    \code{"interact"})
    \item
    a fitted interpoint interaction (object of class
    \code{"fii"})
    \item
    a point process model for simulation (object of class
    \code{"rmhmodel"}), usually obtained from \code{\link[spatstat.random]{rmhmodel}}.
  }
}
\section{Other types of models}{
  Methods for \code{reach} are also defined for
  point process models of class \code{"kppm"} and \code{"dppm"}.
  Their technical definition is different from this one.
  See \code{\link[spatstat.model]{reach.kppm}} and \code{\link[spatstat.model]{reach.dppm}}.
}
\seealso{
    \code{\link[spatstat.model]{reach.ppm}}

     \code{\link{rmhmodel}}

    See \code{\link[spatstat.model]{reach.kppm}} and \code{\link[spatstat.model]{reach.dppm}}
    for other types of point process models.
}
\examples{
    reach(rmhmodel(cif='poisson', par=list(beta=100)))
    # returns 0

    reach(rmhmodel(cif='strauss', par=list(beta=100, gamma=0.1, r=7)))
    # returns 7

    reach(rmhmodel(cif='sftcr', par=list(beta=100, sigma=1, kappa=0.7)))
    # returns Inf
    
    reach(rmhmodel(cif='multihard',
                   par=list(beta=c(10,10), hradii=matrix(c(1,3,3,1),2,2))))
    # returns 3
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{models}
