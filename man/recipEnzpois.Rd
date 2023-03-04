\name{recipEnzpois}
\alias{recipEnzpois}
\title{
  First Reciprocal Moment of the Truncated Poisson Distribution
}
\description{
  Computes the first reciprocal moment (first negative moment)
  of the truncated Poisson distribution
  (the Poisson distribution conditioned to have a nonzero value).
}
\usage{
  recipEnzpois(mu, exact=TRUE)
}
\arguments{
  \item{mu}{
    The mean of the original Poisson distribution.
    A single positive numeric value, or a vector of positive numbers.
  }
  \item{exact}{
    Logical value specifying whether to use the exact analytic formula
    if possible.
  }
}
\details{
  This function calculates the expected value of \eqn{1/N}
  given \eqn{N > 0}, where \eqn{N} is a Poisson random variable
  with mean \eqn{\mu}.

  If the library \pkg{gsl} is loaded, and if \code{exact=TRUE} (the
  default), then the calculation uses
  the exact analytic formula
  \deqn{
    \nu = \frac{e^{-\mu}}{1- e^{-\mu}}
    \left( Ei(\mu) - \log \mu - \gamma \right)
  }{
    nu = (exp(-mu)/(1 - exp(-mu))) (Ei(mu) - log(mu) - gamma)
  }
  where \eqn{\nu} is the desired reciprocal moment, and
  \deqn{
    Ei(x) = \int_{-\infty}^x t e^{-t} dt
  }{
    Ei(x) = int_[-infty]^x t * exp(-t) dt
  }
  is the first exponential integral, and
  \eqn{\gamma \approx 0.577}{\gamma ~= 0.577}
  is the Euler-Mascheroni constant.

  If \pkg{gsl} is not loaded, or if \code{exact=FALSE} is specified,
  the value is computed approximately (and more slowly)
  by summing over the possible values of \eqn{N} up to a finite limit.
}
\value{
  A single numerical value or a numeric vector.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{rpoisnonzero}}
}
\examples{
  if(require(gsl)) {
    v <- recipEnzpois(10)
    print(v)
  }
  recipEnzpois(10, exact=FALSE)
}
\keyword{math}
\keyword{distribution}
