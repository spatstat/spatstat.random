\name{rpoistrunc}
\alias{rpoisnonzero}
\alias{rpoistrunc}
\title{
  Random Values from the Truncated Poisson Distribution
}
\description{
  Generate realisations of a Poisson random variable
  which are truncated, that is, conditioned to be nonzero, or conditioned
  to be at least a given number.
}
\usage{
rpoisnonzero(n, lambda, method=c("harding", "transform"), implem=c("R", "C"))
rpoistrunc(n, lambda, minimum = 1, method=c("harding", "transform"), implem=c("R", "C"))
}
\arguments{
  \item{n}{Number of random values to be generated.}
  \item{lambda}{
    Mean value of the un-truncated Poisson distribution.
    A nonnegative number, or vector of nonnegative numbers.
  }
  \item{minimum}{
    Minimum permitted value for the random variables.
    A nonnegative integer, or vector of nonnegative integers.
  }
  \item{method}{
    Character string (partially matched) specifying the simulation
    algorithm to be used. See Details.
  }
  \item{implem}{
    Character string specifying whether to use the implementation
    in interpreted R code (\code{implem="R"}, the default)
    or C code (\code{implem="C"}).
  }
}
\details{
  \code{rpoisnonzero}
  generates realisations of the
  Poisson distribution with mean \code{lambda}
  conditioned on the event that the values are not equal to zero.

  \code{rpoistrunc} 
  generates realisations of the
  Poisson distribution with mean \code{lambda}
  conditioned on the event that the values are greater than
  or equal to \code{minimum}.
  The default \code{minimum=1} is equivalent to
  generating zero-truncated Poisson random variables
  using \code{rpoisnonzero}.
  The value \code{minimum=0} is equivalent to
  generating un-truncated Poisson random variables
  using \code{\link[stats]{rpois}}.

  The arguments \code{lambda} and \code{minimum} can be
  vectors of length \code{n}, specifying different means
  for the un-truncated Poisson distribution, and different
  minimum values, for each of the \code{n} random output values.

  If \code{method="transform"} the simulated values are generated
  by transforming a uniform random variable using the quantile function
  of the Poisson distribution.
  If \code{method="harding"} (the default) the simulated values are
  generated using an algorithm proposed by E.F. Harding which exploits
  properties of the Poisson point process. The Harding algorithm seems
  to be faster.
}
\value{
  An integer vector of length \code{n}.
}
\author{
  \adrian, after ideas of Ted Harding and Peter Dalgaard.
}
\references{
  E.F. Harding (2005) Email to R-help email group, 01 May 2005.

  P. Dalgaard (2005) Email to R-help email group, 01 May 2005.
}
\seealso{
  \code{\link[stats]{rpois}} for Poisson random variables.

  \code{\link{recipEnzpois}} for the reciprocal moment of \code{rpoisnonzero}.
}
\examples{
  rpoisnonzero(10, 0.8)

  rpoistrunc(10, 1, 2)
}
\keyword{datagen}
