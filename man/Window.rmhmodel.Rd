\name{Window.rmhmodel}
\alias{Window.rmhmodel}
\title{Extract Window of Spatial Object}
\description{
  Given a spatial object (such as a point pattern or pixel image)
  in two dimensions, these functions extract the window
  in which the object is defined.
}
\usage{


 \method{Window}{rmhmodel}(X, \dots)


}
\arguments{
  \item{X}{A spatial object.}
  \item{\dots}{Ignored.}
}
\value{
  An object of class \code{"owin"} (see \code{\link[spatstat.geom]{owin.object}})
  specifying an observation window.
}
\details{
  These are methods for the generic function \code{\link[spatstat.geom]{Window}}
  which extract the spatial window in which the object \code{X}
  is defined. 
}
\seealso{
  \code{\link[spatstat.geom]{Window}},
  \code{\link[spatstat.geom]{Window.ppp}},
  \code{\link[spatstat.geom]{Window.psp}}.
  
  \code{\link[spatstat.geom]{owin.object}}
}
\examples{
   A <- rmhmodel(cif='poisson', par=list(beta=10), w=square(2))
   Window(A)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
 