\name{rdiffuse}
\alias{rdiffuse}
\alias{rdiffuse.ppp}
\title{
  Perturb the Points in a Point Pattern According to a Diffusion Process
}
\description{
  Given a spatial point pattern inside a window, allow each point of the pattern
  to undergo random spatial diffusion inside the window
  for a specified amount of time,
  and return the final position of each point.
}
\usage{
rdiffuse(X, sigma, \dots)

\method{rdiffuse}{ppp}(X, sigma, \dots, nsim=1, drop=TRUE,
    connect = 8, method = c("C", "interpreted"), unround = TRUE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"})
    specifying the initial position of each point.
  }
  \item{sigma}{
    Equivalent standard deviation of the final position of each point,
    if the window were unbounded. A single positive number,
    or a pixel image.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}}
    controlling the spatial resolution.
  }
  \item{nsim}{
    Number of simulated realisations to be generated.
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
  \item{connect}{
    Pixel grid connectivity (4 or 8).
  }
  \item{method}{
    For testing purposes only.
    Character string (partially matched) specifying whether to
    use a \proglang{C} language implementation or an interpreted \R
    implementation. 
  }
  \item{unround}{
    Logical value specifying whether the final positions
    of the points should be altered slightly 
    by adding a small random error to the coordinates
    using \code{\link[spatstat.random]{rUnround}}. This reverses the effect of
    the discretisation.
    If \code{unround=FALSE}, all the final positions are exactly at the
    centre of a pixel.    
  }
}
\details{
  The spatial locations of the points of \code{X} are first discretised
  onto a pixel grid, with resolution specified by the arguments
  \code{\dots}.
  Each point then executes a random walk on the pixel grid,
  independently of other points. The final location of each point
  is returned.
}
\value{
  A point pattern or a list of point patterns.
  Each point pattern has the same window as \code{X} and the same number of
  points as \code{X}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.explore]{densityHeat.ppp}}
}
\examples{
  plot(solist(original=cells, diffused=rdiffuse(cells, 0.05)))
}
\keyword{datagen}
\concept{Diffusion}

