\name{methods.clusterprocess}
\alias{methods.clusterprocess} % DoNotExport
\alias{predict.clusterprocess}
\alias{intensity.clusterprocess}
\alias{print.clusterprocess}
\alias{clusterradius.clusterprocess}
\alias{reach.clusterprocess}
\alias{simulate.clusterprocess}
\title{
  Methods for Cluster Models
}
\description{
  Methods for the class \code{"clusterprocess"} of cluster point process
  models with specified values of the parameters.
}
\usage{
 \method{intensity}{clusterprocess}(X, \dots)

 \method{predict}{clusterprocess}(object, \dots,
                  locations, type = "intensity", ngrid = NULL)

 \method{print}{clusterprocess}(x, \dots)

 \method{clusterradius}{clusterprocess}(model,\dots,thresh=NULL, precision=FALSE)

 \method{reach}{clusterprocess}(x, \dots, epsilon)

 \method{simulate}{clusterprocess}(object, nsim=1, \dots, win=unit.square(), window=win)
}
\arguments{
  \item{model,object,x,X}{
    Object of class \code{"clusterprocess"}.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{locations}{
    Locations where prediction should be performed.
    A window or a point pattern.
  }
  \item{type}{
    Currently must equal \code{"intensity"}.
  }
  \item{ngrid}{
    Pixel grid dimensions for prediction, if \code{locations} is
    a rectangle or polygon.
  }
  \item{thresh,epsilon}{Tolerance thresholds}
  \item{precision}{
    Logical value stipulating whether the precision should also be returned.
  }
  \item{win,window}{
    Window (object of class \code{"owin"}) in which the simulated
    pattern should be generated.
  }
  \item{nsim}{
    Number of simulated patterns to be generated.
  }
}
\details{
  An object of class \code{"clusterprocess"} represents a Neyman-Scott-Cox
  cluster process model with specified parameter values.

  This page documents methods for printing, simulating and predicting
  such models, supplied by the package \pkg{spatstat.random}.
  Other methods are supplied in other packages.
}
\value{
  Same as for other methods.
}
\author{
  \adrian
}
\seealso{
  \code{\link{clusterprocess}}
  %% \code{\link[spatstat.model]{Kmodel.clusterprocess}} WHEN AVAILABLE
}
\examples{
  m <- clusterprocess("Thomas", kappa=10, mu=5, scale=0.1)
  m2 <- clusterprocess("VarGamma", kappa=10, mu=10, scale=0.1, nu=0.7)
  m
  m2
  intensity(m)
  reach(m)
  clusterradius(m)
  Z <- predict(m, locations=square(2))
  X <- simulate(m)
}
\keyword{spatial}
\keyword{models}
