\name{domain.rmhmodel}
\alias{domain.rmhmodel}
\title{
  Extract the Domain of any Spatial Object
}
\description{
  Given a spatial object such as a point pattern, in any number of dimensions,
  this function extracts the spatial domain in which the object is defined.
}
\usage{




 \method{domain}{rmhmodel}(X, \dots)


}
\arguments{
  \item{X}{
    A spatial object such as a point pattern (in any number
    of dimensions), line segment pattern or pixel image.
  }
  \item{\dots}{
    Extra arguments. They are ignored by all the methods listed here.
  }
}
\details{
  The function \code{\link[spatstat.geom]{domain}} is generic.

  For a spatial object \code{X} in any number of dimensions, 
  \code{domain(X)} extracts the spatial domain in which \code{X} is
  defined.

  For a two-dimensional object \code{X}, typically \code{domain(X)}
  is the same as \code{Window(X)}.

  Exceptions occur for methods related to linear networks.

}
\value{
  A spatial object representing the domain of \code{X}.
  Typically a window (object of class \code{"owin"}),
  a three-dimensional box (\code{"box3"}), a multidimensional
  box (\code{"boxx"}) or a linear network (\code{"linnet"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.geom]{domain}},
  \code{\link[spatstat.geom]{domain.quadratcount}},
  \code{\link[spatstat.model]{domain.ppm}},
  \code{\link[spatstat.explore]{domain.quadrattest}},
  \code{\link[spatstat.linnet]{domain.lpp}}.
  \code{\link[spatstat.geom]{Window}},
  \code{\link[spatstat.geom]{Frame}}.
}
\examples{
  domain(rmhmodel(cif='poisson', par=list(beta=1), w=square(2)))
}
\keyword{spatial}
\keyword{manip}
