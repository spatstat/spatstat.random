\name{dpakes}
\alias{dpakes}
\alias{ppakes}
\alias{qpakes}
\alias{rpakes}
\title{
  Pakes distribution
}
\description{
  Probability density, cumulative distribution function,
  quantile function, and random generation for the Pakes
  distribution.
}
\usage{
dpakes(x, zeta)
ppakes(q, zeta)
qpakes(p, zeta)
rpakes(n, zeta)
}
\arguments{
  \item{x,q}{
    Numeric vector of quantiles.
  }
  \item{p}{Numeric vector of probabilities}
  \item{n}{Number of observations.}
  \item{zeta}{
    Mean of distribution. A single, non-negative, numeric value.
  }
}
\details{
  These functions concern the probability distribution of the random
  variable
  \deqn{
    X = \sum_{n=1}^\infty \prod_{j=1}^n U_j^{1/\zeta}
  }{
    X = sum_[n=1]^[infty] prod_[j=1]^[n] U[j]^(1/zeta)
  }
  where \eqn{U_1, U_2, \ldots}{U[1], U[2], ...} are independent
  random variables uniformly distributed on \eqn{[0,1]}
  and \eqn{\zeta}{zeta} is a parameter.

  This distribution arises in many contexts. For example, 
  for a homogeneous Poisson point process in two-dimensional space
  with intensity \eqn{\lambda}{lambda},
  the standard Gaussian kernel estimator of intensity
  with bandwidth \eqn{\sigma}, evaluated at any fixed location \eqn{u},
  has the same distribution as \eqn{(\lambda/\zeta) X}{(lambda/zeta) * X}
  where
  \eqn{\zeta = 2 \pi \lambda\sigma^2}{zeta = 2 * pi * lambda * sigma^2}.
  
  Following the usual convention,
  \code{dpakes} computes the probability density,
  \code{ppakes} the cumulative distribution function,
  and \code{qpakes} the quantile function,
  and \code{rpakes} generates random variates with this distribution.

  The computation is based on a recursive integral equation
  for the cumulative distribution function, due to Professor Tony Pakes,
  presented in Baddeley, Moller and Pakes (2008).
  The solution uses the fact that the random variable
  satisfies the distributional equivalence
  \deqn{
    X \equiv U^{1/\zeta} (1 + X)
  }{
    X == U^(1/zeta) * (1+X)
  }
  where \eqn{U} is uniformly distributed on \eqn{[0,1]} and independent
  of \eqn{X}.
}
\value{
  A numeric vector. 
}
\references{
  Baddeley, A.,  Moller, J. and Pakes, A.G. (2008)
  Properties of residuals for spatial point processes,
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{60}, 627--649.
}
\author{
  Adrian Baddeley.
}
\examples{
 curve(dpakes(x, 1.5), to=4)
 rpakes(3, 1.5)
}
\keyword{distribution}
