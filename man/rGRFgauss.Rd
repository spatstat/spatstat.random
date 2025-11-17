\name{rGRFgauss}
\alias{rGRFgauss}
\alias{rGRFexpo}
\alias{rGRFstable}
\alias{rGRFgencauchy}
\alias{rGRFmatern}
\title{
  Simulate a Gaussian Random Field
}
\description{
  Generate a realisation of a Gaussian Random Field
  with given first and second moments.
}
\usage{
rGRFgauss(W = owin(), mu = 0, var = 1, scale,
          \dots, nsim = 1, drop = TRUE)
rGRFexpo(W = owin(), mu = 0, var = 1, scale,
          \dots, nsim = 1, drop = TRUE)
rGRFstable(W = owin(), mu = 0, var = 1, scale, alpha,
           \dots, nsim = 1, drop = TRUE)
rGRFgencauchy(W = owin(), mu = 0, var = 1, scale, alpha, beta, 
           \dots, nsim = 1, drop = TRUE)
rGRFmatern(W = owin(), mu = 0, var = 1, scale, nu, 
           \dots, nsim = 1, drop = TRUE)
}
\arguments{
  \item{W}{
    Window (object of class \code{"owin"}) in which to generate the
    simulated random field.
  }
  \item{mu}{
    The mean of the random field. A single numeric value,
    or a \code{function(x, y, ...)},
    or a pixel image (object of class \code{"im"}).
  }
  \item{var}{
    Variance of the random field. A single positive number.
  }
  \item{scale}{
    Spatial scale parameter \eqn{h}. A single positive number.
    See Details.
  }
  \item{alpha,beta,nu}{
    Additional parameters for specific models. See Details.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} to determine
    the pixel resolution.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  These functions generate simulated realisations
  of a Gaussian random field.

  The mean \eqn{E[Z(u)]} of the Gaussian random field value \eqn{Z(u)}
  at any location \eqn{u} is specified by the argument \code{mu},
  which may be a constant, a \code{function(x,y,...)}, or a pixel image.

  The variance \eqn{V[Z(u)]} of the Gaussian random field value
  is specified by the argument \code{var}, which should be a single
  positive numerical value.

  The correlation \eqn{C(u - v) = C(Z(u), Z(v))} between the values
  at two locations \eqn{u} and \eqn{v} depends on the
  distance \eqn{r = \| u-v\|}{r = |u-v|} as follows:
  
  \describe{
    \item{\code{rGRFexpo}}{
      the exponential covariance function
      \deqn{C(r) = \sigma^2 \exp(-r/h)}{C(r) = sigma^2 * exp(-r/h)}
      where \eqn{\sigma^2} is the variance parameter \code{var},
      and \eqn{h} is the scale parameter \code{scale}.
    }
    \item{\code{rGRFgauss}}{
      the Gaussian covariance function
      \deqn{C(r) = \sigma^2 \exp(-(r/h)^2)}{C(r) = sigma^2 * exp(-(r/h)^2)}
      where \eqn{\sigma^2} is the variance parameter \code{var},
      and \eqn{h} is the scale parameter \code{scale}.
    }
    \item{\code{rGRFstable}}{
      the stable covariance function
      \deqn{
	C(r) = \sigma^2 \exp(-(r/h)^\alpha)
      }{
	C(r) = sigma^2 * exp(-(r/h)^alpha)
      }
      where \eqn{\sigma^2} is the variance parameter \code{var},
      \eqn{h} is the scale parameter \code{scale},
      and \eqn{\alpha} is the shape parameter \code{alpha}.
    }
    \item{\code{rGRFgencauchy}}{
      the generalised Cauchy covariance function
      \deqn{
	C(r) = \sigma^2 (1 + (x/h)^\alpha)^{-\beta/\alpha}
      }{
	C(r) = sigma^2 * (1 + (x/h)^\alpha)^(-\beta/\alpha)
      }
      where \eqn{\sigma^2} is the variance parameter \code{var},
      \eqn{h} is the scale parameter \code{scale},
      and \eqn{\alpha} and \eqn{\beta} are the shape parameters
      \code{alpha} and \code{beta}.
    }
    \item{\code{rGRFmatern}}{
      the Whittle-\Matern covariance function
      \deqn{
	C(r) = \sigma^2 \frac{1}{2^{\nu-1} \Gamma(\nu)}
	(\sqrt{2 \nu} \, r/h)^\nu K_\nu(\sqrt{2\nu}\, r/h)
      }{
	C(r) = \sigma^2 * 2^(1-\nu) * \Gamma(\nu)^(-1) *
	(sqrt(2 *\nu) * r/h)^\nu * K[\nu](sqrt(2 * nu) * r/h)
      }
      where \eqn{\sigma^2} is the variance parameter \code{var},
      \eqn{h} is the scale parameter \code{scale},
      and \eqn{\nu} is the shape parameter \code{nu}.
    }
  }
  The algorithm generates \code{nsim} simulated realisations
  of the random field using the circulant embedding technique
  (Davies and Bryant, 2013).
}
\value{
  If \code{nsim=1} and \code{drop=TRUE}, a pixel image
  (object of class \code{"im"}).
  Otherwise, a list of pixel images.
}
\references{
  Davies, T.M. and Bryant, D. (2013)
  On Circulant Embedding for Gaussian Random Fields in R.
  \emph{Journal of Statistical Software} 55, issue 9
  \doi{10.18637/jss.v055.i09}
}
\author{
  \tilman and David Bryant.
  Modified by \adrian.
}
\seealso{
  \code{\link[spatstat.random]{rLGCP}}
}
\examples{
  plot(rGRFgauss(scale=0.1))
}
\keyword{spatial}
\keyword{datagen}
