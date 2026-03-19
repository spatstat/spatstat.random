\name{rpoint3}
\alias{rpoint3}
\title{
  Generate N Random Points in 3 Dimensions
}
\description{
  Generate a random point pattern in three-dimensional space
  containing \eqn{n} independent, identically distributed random points
  with a specified distribution.
}
\usage{
rpoint3(n, f, fmax = NULL, domain = box3(), \dots,
  nsim = 1, drop = TRUE, ex = NULL, giveup = 1000)
}
\arguments{
  \item{n}{
    Number of points to generate.
  }
  \item{f}{
    The probability density of the points,
    possibly un-normalised.
    Either a constant, or a function \code{f(x,y,z, \dots)}.
  }
  \item{fmax}{
    An upper bound on the values of \code{f}.
  }
  \item{domain}{
    Three-dimensional box (object of class \code{"box3"})
    in which to generate the point pattern.
  }
  \item{\dots}{
    Additional arguments passed to the function \code{f}.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
  \item{ex}{
    Optional. Example point pattern in three dimensions
    (object of class \code{"pp3"}). The number of points \code{n}
    and the spatial region \code{domain} will be taken from \code{ex}.
  }
  \item{giveup}{
    Integer. Maximum number of attempts permitted in the rejection method.
  }
}
\details{
  This function generates \code{n} independent, identically distributed
  random points with common probability density proportional to
  \code{f} inside the specified three-dimensional region \code{domain}.

  The argument \code{f} may be
  \describe{
    \item{a numerical constant:}{
      uniformly distributed random points will
      be generated.
    }
    \item{a function:}{random points will be generated
      in the \code{domain} with probability density proportional
      to \code{f(x, y, z, \dots)} where \code{x}, \code{y}, \code{z}
      are the cartesian
      coordinates. The function \code{f} must accept 
      \emph{vectors} of coordinates \code{x,y,z} and return the corresponding
      vector of function values. Additional arguments \code{\dots} of any kind
      may be passed to the function.
    }
  }
  The algorithm is as follows:
  \itemize{
    \item
    If \code{f} is a constant, we invoke \code{\link{runifpoint3}}.
    \item
    If \code{f} is a function, then we use the rejection method.
    Proposal points are generated from the uniform distribution.
    A proposal point \eqn{(x,y,z)} is accepted with probability
    \code{f(x,y,z,\dots)/fmax} and otherwise rejected.
    The algorithm continues until \code{n} points have been
    accepted. It gives up after \code{giveup * n} proposals
    if there are still fewer than \code{n} points.
  }
}
\value{
  A three-dimensional point pattern (an object of class \code{"pp3"})
  if \code{nsim=1} and \code{drop=TRUE},
  otherwise a list of three-dimensional point patterns.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{runifpoint3}}
}
\examples{
  rpoint3(10, 1)
  rpoint3(10, f=function(x,y,z) {sqrt(x)}, fmax=1)
}
\keyword{spatial}
\keyword{datagen}
\concept{Three-dimensional}
