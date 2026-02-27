\name{clusterprocess}
\alias{clusterprocess}
\title{
  Cluster Point Process Model
}
\description{
  Creates an object representing a cluster point process model
  with the specified parameters.
  Typically used for simulations or calculations about such a model.
}
\usage{
clusterprocess(name = "Thomas", \dots, mu, kappa, scale)
}
\arguments{
  \item{name}{
    Name of the cluster process. One of
    \code{"Thomas"},
    \code{"MatClust"},
    \code{"VarGamma"} or
    \code{"Cauchy"}.
  }
  \item{\dots}{
    Other arguments needed for the model.
  }
  \item{mu}{
    Mean cluster size. A single number, or a pixel image.
  }
  \item{kappa}{
    Parent intensity. A single number.
  }
  \item{scale}{
    Cluster scale parameter of the model.
  }
}
\details{
  This function creates an object representing a Neyman-Scott-Cox
  cluster process model with the specified parameter values.
  The object belongs to the class \code{"clusterprocess"}.

  This is different from the model-fitting function
  \code{\link[spatstat.model]{kppm}} which fits such a model to a point pattern dataset.

  There are methods for simulating a \code{"clusterprocess"},
  calculating its moments, and other purposes, listed under
  \emph{See Also}.
}
\value{
  Object of class \code{"clusterprocess"}.
}
\author{
  \adrian
}
\seealso{
  Methods for 
  \code{simulate},
  \code{predict},
  \code{intensity},
  \code{print},
  \code{clusterradius} and 
  \code{reach}
  are supplied in package \pkg{spatstat.random}
  and documented under \code{\link{methods.clusterprocess}}.

  The functions
  \code{\link[spatstat.model]{Kmodel}},
  \code{\link[spatstat.model]{pcfmodel}},
  \code{\link[spatstat.model]{psib}},
  \code{\link[spatstat.model]{persist}},
  \code{\link[spatstat.model]{clusterstrength}} and
  \code{\link[spatstat.model]{varcount}}
  in package \pkg{spatstat.model}
  can also be applied to a \code{"clusterprocess"}.
}
\examples{
  m <- clusterprocess("Thomas", kappa=10, mu=5, scale=0.1)
  simulate(m, win=square(2))
}
\keyword{spatial}
\keyword{models}

