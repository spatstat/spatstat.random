\name{indefinteg}
\alias{indefinteg}
\title{
  Indefinite Integral
}
\description{
  Computes the indefinite integral of the given function.
}
\usage{
  indefinteg(f, x, \dots,
             method=c("trapezoid", "quadrature"),
             start=min(x), nfine=8192) 
}
\arguments{
  \item{f}{
    an \R function taking a numeric first argument and returning a
    numeric vector of the same length.
  }
  \item{x}{
    Vector of values of the argument for which the indefinite integral
    should be evaluated. 
  }
  \item{\dots}{
    additional arguments to be passed to \code{f}.
  }
  \item{method}{
    String (partially matched) specifying how to compute the integrals.
  }
  \item{start}{
    Lower limit of integration. A single number.
  }
  \item{nfine}{
    Number of sub-intervals to use for computation
    if \code{method='trapezoid'}.
  }
}
\details{
  The indefinite integral of the given function \code{f}
  is computed numerically at each of the desired values \code{x}.
  The lower limit of integration is taken to be \code{min(x)}.

  If \code{method='trapezoid'} (the default),
  the integrals are computed rapidly using the trapezoid rule.
  If \code{method='quadrature'} the integrals are computed
  accurately but much more slowly, using the numerical quadrature routine
  \code{\link[stats]{integrate}}.
  
  If \code{method='trapezoid'} 
  the function \code{f} is first evaluated 
  on a finer grid of values of the function argument.
  The fine grid contains at least \code{nfine} sample points.
  The values of the indefinite integral on the fine grid
  are computed using the trapezoidal approximation.
  Finally the values of the indefinite integral are extracted at
  the desired argument values \code{x}.
}
\value{
  Numeric vector of the same length as \code{x}. 
}
\author{
  \adrian.
}
\seealso{
  \code{\link[stats]{integrate}}
}
\examples{
  curve(indefinteg(sin, x), to=pi)
}
\keyword{math}
