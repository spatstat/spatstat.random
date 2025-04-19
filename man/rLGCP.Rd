\name{rLGCP}
\alias{rLGCP}
\title{Simulate Log-Gaussian Cox Process}
\description{
  Generate a random point pattern, a realisation of the
  log-Gaussian Cox process.
}
\usage{
 rLGCP(model=c("exponential", "gauss", "stable", "gencauchy", "matern"),
       mu = 0, param = NULL,
       \dots,
       win=NULL, saveLambda=TRUE, nsim=1, drop=TRUE,
       n.cond=NULL, w.cond=NULL)
}
\arguments{
  \item{model}{
    character string (partially matched) giving the name of a
    covariance model for the Gaussian random field. 
  }
  \item{mu}{
    mean function of the Gaussian random field. Either a
    single number, a \code{function(x,y, ...)} or a pixel
    image (object of class \code{"im"}).
  }
  \item{param}{
    List of parameters for the covariance.
    Standard arguments are \code{var} and \code{scale}.
  }
  \item{\dots}{
    Additional parameters for the covariance,
    or arguments passed to \code{\link[spatstat.geom]{as.mask}} to determine
    the pixel resolution.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}.
  }
  \item{saveLambda}{
    Logical. If \code{TRUE} (the default) then the
    simulated random intensity will also be saved,
    and returns as an attribute of the point pattern.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
  \item{n.cond}{
    Optional. Integer specifying a fixed number of points.
    See the section on \emph{Conditional Simulation}.
  }
  \item{w.cond}{
    Optional. Conditioning region. A window (object of class \code{"owin"})
    specifying the region which must contain exactly \code{n.cond} points.
    See the section on \emph{Conditional Simulation}.
  }
}
\value{
  A point pattern (object of class \code{"ppp"})
  or a list of point patterns.

  Additionally, the simulated intensity function for each point pattern is
  returned as an attribute \code{"Lambda"} of the point pattern,
  if \code{saveLambda=TRUE}.
}
\details{
  This function generates a realisation of a log-Gaussian Cox
  process (LGCP). This is a Cox point process in which
  the logarithm of the random intensity is a Gaussian random
  field with mean function \eqn{\mu} and covariance function
  \eqn{c(r)}. Conditional on the random intensity, the point process
  is a Poisson process with this intensity.

  The string \code{model} specifies the covariance 
  function of the Gaussian random field, and the parameters
  of the covariance are determined by \code{param} and \code{\dots}.

  All models recognise the parameters \code{var}
  for the variance at distance zero, and \code{scale} for the scale
  parameter. Some models require additional parameters which are listed
  below. 
  
  The available models are as follows:
  \describe{
    \item{\code{model="exponential"}:}{
      the exponential covariance function
      \deqn{C(r) = \sigma^2 \exp(-r/h)}{C(r) = sigma^2 * exp(-r/h)}
      where \eqn{\sigma^2} is the variance parameter \code{var},
      and \eqn{h} is the scale parameter \code{scale}.
    }
    \item{\code{model="gauss"}:}{
      the Gaussian covariance function
      \deqn{C(r) = \sigma^2 \exp(-(r/h)^2)}{C(r) = sigma^2 * exp(-(r/h)^2)}
      where \eqn{\sigma^2} is the variance parameter \code{var},
      and \eqn{h} is the scale parameter \code{scale}.
    }
    \item{\code{model="stable"}:}{
      the stable covariance function
      \deqn{
	C(r) = \sigma^2 \exp(-(r/h)^\alpha)
      }{
	C(r) = sigma^2 * exp(-(r/h)^alpha)
      }
      where \eqn{\sigma^2} is the variance parameter \code{var},
      \eqn{h} is the scale parameter \code{scale},
      and \eqn{\alpha} is the shape parameter \code{alpha}.
      The parameter \code{alpha} must be given, either as a stand-alone
      argument, or as an entry in the list \code{param}.
    }
    \item{\code{model="gencauchy"}:}{
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
      The parameters \code{alpha} and \code{beta}
      must be given, either as stand-alone arguments, or as entries
      in the list \code{param}.
    }
    \item{\code{model="matern"}:}{
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
      The parameter \code{nu} must be given, either as a stand-alone
      argument, or as an entry in the list \code{param}.
    }
  }
  The algorithm uses the circulant embedding technique to 
  generate values of a Gaussian random field,
  with the specified mean function \code{mu}
  and the covariance specified by the arguments \code{model} and
  \code{param}, on the points of a regular grid. The exponential
  of this random field is taken as the intensity of a Poisson point
  process, and a realisation of the Poisson process is then generated by the 
  function \code{\link[spatstat.random]{rpoispp}} in the \pkg{spatstat.random} package.
  
  If the simulation window \code{win} is missing or \code{NULL},
  then it defaults to 
  \code{Window(mu)} if \code{mu} is a pixel image,
  and it defaults to the unit square otherwise.
  
  The LGCP model can be fitted to data using \code{\link[spatstat.model]{kppm}}.
}
\section{Conditional Simulation}{
  If \code{n.cond} is specified, it should be a single integer.
  Simulation will be conditional on the event 
  that the pattern contains exactly \code{n.cond} points
  (or contains exactly \code{n.cond} points inside
  the region \code{w.cond} if it is given).
 
  Conditional simulation uses the rejection algorithm described
  in Section 6.2 of Moller, Syversveen and Waagepetersen (1998).
  There is a maximum number of proposals which will be attempted.
  Consequently the return value may contain fewer
  than \code{nsim} point patterns.
}
\section{Warning: new implementation}{
  The simulation algorithm for \code{rLGCP} has been completely re-written
  in \pkg{spatstat.random} version \code{3.2-0} to avoid depending on
  the package \pkg{RandomFields} which is now defunct (and is sadly missed).

  It is no longer possible to replicate results that were obtained using
  \code{rLGCP} in previous versions of \pkg{spatstat.random}.

  The current code is a new implementation and should be considered
  vulnerable to new bugs. 
}
\seealso{
\code{\link[spatstat.random]{rpoispp}},
\code{\link[spatstat.random]{rMatClust}},
\code{\link[spatstat.random]{rGaussPoisson}},
\code{\link[spatstat.random]{rNeymanScott}}.

For fitting the model, see \code{\link[spatstat.model]{kppm}},
\code{\link[spatstat.model]{lgcp.estK}}.
}
\references{
  \Moller, J., Syversveen, A. and Waagepetersen, R. (1998)
  Log Gaussian Cox Processes.
  \emph{Scandinavian Journal of Statistics} \bold{25}, 451--482.
}
\examples{
    online <- interactive()

    # homogeneous LGCP with exponential covariance function
    X <- rLGCP("exp", 3, var=0.2, scale=.1)

    # inhomogeneous LGCP with Gaussian covariance function
    m <- as.im(function(x, y){5 - 1.5 * (x - 0.5)^2 + 2 * (y - 0.5)^2}, W=owin())
    X <- rLGCP("gauss", m, var=0.15, scale =0.1)

    if(online) {
     plot(attr(X, "Lambda"))
     points(X)
    }

    # inhomogeneous LGCP with Matern covariance function
    X <- rLGCP("matern", function(x, y){ 1 - 0.4 * x},
               var=2, scale=0.7, nu=0.5,
               win = owin(c(0, 10), c(0, 10)))
    if(online) plot(X)
}
\author{Abdollah Jalilian and Rasmus Waagepetersen.
  Modified by \spatstatAuthors.
}
\keyword{spatial}
\keyword{datagen}

