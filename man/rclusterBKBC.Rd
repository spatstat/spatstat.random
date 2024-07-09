\name{rclusterBKBC}
\alias{rclusterBKBC}
\title{
  Simulate Cluster Process using Brix-Kendall Algorithm or Modifications
}
\description{
  Generates simulated realisations of a stationary Neyman-Scott
  cluster point process, using the Brix-Kendall (2002) algorithm or various
  modifications proposed by Baddeley and Chang (2023).
  For advanced research use.
}
\usage{
rclusterBKBC(clusters="Thomas",
   kappa, mu, scale,
   \dots,
   W = unit.square(),
   nsim = 1, drop = TRUE,
   best = FALSE,
   external = c("BK", "superBK", "border"),
   internal = c("dominating", "naive"),
   inflate = 1,
   psmall = 1e-04,
   use.inverse=TRUE,
   use.special=TRUE,
   integralmethod=c("quadrature", "trapezoid"),
   verbose = TRUE, warn=TRUE)
}
\arguments{
  \item{clusters}{
    Character string (partially matched) specifying the cluster process.
    Current options include \code{"Thomas"}, \code{"MatClust"},
    \code{"Cauchy"} and \code{"VarGamma"}.
  }
  \item{kappa}{
    Intensity of the parent process. A nonnegative number.
  }
  \item{mu}{
    Mean number of offspring per parent. A nonnegative number.
  }
  \item{scale}{
    Cluster scale. Interpretation depends on the model.
  }
  \item{\dots}{
    Additional arguments controlling the shape of the cluster kernel, if any.
  }
  \item{W}{
    Window in which the simulation should be generated.
    An object of class \code{"owin"}.
  }
  \item{nsim}{
    The number of simulated point patterns to be generated. A positive integer.
  }
  \item{drop}{
    Logical value. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
  \item{best}{
    Logical value. If \code{best=TRUE}, the code will choose the fastest
    algorithm. If \code{best=FALSE} (the default), the algorithm will be
    specified by the other arguments \code{external} and
    \code{internal}.
    See Details.
  }
  \item{external}{
    Algorithm to be used to generate parent points which lie outside the
    bounding window. See Details. 
  }
  \item{internal}{
    Algorithm to be used to generate parent points which lie inside the
    bounding window. See Details.
  }
  \item{inflate}{
    Numerical value determining the position of the bounding window.
    See Details.
  }
  \item{psmall}{
    Threshold of small probability for use in the algorithm.
  }
  \item{use.inverse}{
    Logical value specifying whether to compute the inverse function
    analytically, if possible (\code{use.inverse=TRUE}, the default)
    or by numerical root-finding (\code{use.inverse=FALSE}).
    This is mainly for checking validity of code.
  }
  \item{use.special}{
    Logical value specifying whether to use efficient special code
    (if available) to generate the simulations (\code{use.special=TRUE},
    the default) or to use generic code (\code{use.special=FALSE}).
    This is mainly for checking validity of code.
  }
  \item{integralmethod}{
    Character string (partially matched)
    specifying how to perform numerical computation of integrals
    when required. This argument is passed to
    \code{\link[spatstat.univar]{indefinteg}}.
    The default \code{integralmethod="quadrature"} is accurate but
    can be slow. 
    Faster, but possibly less accurate, integration can be performed
    by setting \code{integralmethod="trapezoid"}.
  }
  \item{verbose}{
    Logical value specifying whether to print detailed information about
    the simulation algorithm during execution.
  }
  \item{warn}{
    Logical value specifying whether to issue a warning
    if the number of random proposal points is very large.
  }
}
\details{
  This function is intended for advanced research use.
  It implements the algorithm of Brix and Kendall (2002)
  for generating simulated realisations of a stationary Neyman-Scott
  process, and various modifications of this algorithm proposed
  in Baddeley and Chang (2023).
  It is an alternative to \code{\link[spatstat.random]{rNeymanScott}}.

  The function supports the following models:
  \itemize{
    \item \code{clusters="Thomas"}: the (modified) Thomas cluster process
    which can also be simulated by
    \code{\link[spatstat.random]{rThomas}}.  
    \item \code{clusters="MatClust"}: the \Matern cluster process
    which can also be simulated by
    \code{\link[spatstat.random]{rMatClust}}.  
    \item \code{clusters="Cauchy"}: the Cauchy cluster process
    which can also be simulated by
    \code{\link[spatstat.random]{rCauchy}}.  
    \item \code{clusters="VarGamma"}: the variance-gamma cluster process
    which can also be simulated by
    \code{\link[spatstat.random]{rVarGamma}}.  
    \item 
    any other Poisson cluster process models that may be recognised by 
    \code{\link[spatstat.model]{kppm}}. 
  }

  By default, the code executes the original Brix-Kendall algorithm
  described in Sections 2.3 and 3.1 of Brix and Kendall (2002).

  Modifications of this algorithm, proposed in Baddeley and Chang
  (2023), can be selected using the
  arguments \code{external} and \code{internal}, or \code{best}.

  If \code{best=TRUE}, the code will choose the algorithm
  that would run fastest with the given parameters.
  If \code{best=FALSE} (the default), the choice of algorithm
  is determined by the arguments \code{external} and \code{internal}.

  First the window \code{W} is enclosed in a disc \code{D}
  and Monte Carlo proposal densities are defined with reference to \code{D}
  as described in Brix and Kendall (2002).
  Then \code{D} is inflated by the scale factor \code{inflate}
  to produce a larger disc \code{E} (by default \code{inflate=1}
  implying \code{E=D}). 
  Then the parent points of the clusters are generated, possibly
  using different mechanisms inside and outside \code{E}.

  The argument \code{external} determines the algorithm for generating
  parent points outside \code{E}.
  \itemize{
    \item
    If \code{external="BK"} (the default), proposed parents outside
    \code{E} will be generated from a dominating point process as described in
    Section 3.1 of Brix and Kendall (2002). These points will be thinned
    to obtain the correct intensity of parent points.
    For each accepted parent, offspring points are generated inside \code{D},
    subject to the condition that the parent has at least one offspring inside \code{D}.
    Offspring points are subsequently clipped to the true window \code{W}.
    \item 
    If \code{external="superBK"}, proposed parents will initially be generated
    from a process that dominates the dominating point process
    as described in Baddeley and Chang (2023).
    These proposals will then be thinned to obtain the correct intensity
    of the dominating process, then thinned again to obtain the correct
    intensity of parent points. This procedure reduces computation time
    when \code{scale} is large.
    For each accepted parent, offspring points are generated inside \code{D},
    subject to the condition that the parent has at least one offspring inside \code{D}.
    Offspring points are subsequently clipped to the true window \code{W}.
    \item
    If \code{external="border"} then proposed parents will be generated
    with uniform intensity in a border region
    surrounding the disc \code{D}.
    For each proposed parent, offspring points are generated in the
    entire plane according to the cluster offspring distribution, without
    any restriction.
    Offspring points are subsequently clipped to the true window \code{W}.
    This is the technique currently used in
    \code{\link[spatstat.random]{rNeymanScott}}.
  }
  The argument \code{internal} determines the algorithm for generating
  proposed parent points inside \code{E}.
  \itemize{
    \item
    If \code{internal="dominating"}, parent points in \code{E} are generated
    according to the dominating point process described in 
    Sections 2.3 and 3.1 of Brix and Kendall (2002), and then thinned
    to obtain the correct intensity of parent points. 
    For each accepted parent, offspring points are generated inside \code{D},
    subject to the condition that the parent has at least one offspring inside \code{D}.
    Offspring points are subsequently clipped to the true window \code{W}.
    \item
    If \code{internal="naive"}, parent points in \code{E} are generated with
    uniform intensity inside \code{E} and are not thinned.
    For each proposed parent, offspring points are generated in the
    entire plane according to the cluster offspring distribution, without
    any restriction.
    Offspring points are subsequently clipped to the true window
    \code{W}.
    This is the technique currently used in \code{\link[spatstat.random]{rNeymanScott}}.
  }
  If \code{warn=TRUE}, then a warning will be issued if
  the number of random proposal points (proposed parents and proposed
  offspring) is very large.
  The threshold is \code{\link[spatstat.geom]{spatstat.options}("huge.npoints")}.
  This warning has no consequences,
  but it helps to trap a number of common problems.
}
\value{
  A point pattern, or a list of point patterns.

  If \code{nsim=1} and \code{drop=TRUE}, the result is 
  a point pattern (an object of class \code{"ppp"}).

  Otherwise, the result is a list of \code{nsim} point patterns,
  and also belongs to the class \code{"solist"}.
}
\author{
  \adrian and \yamei.
}
\references{
  \baddchangclustersim

  Brix, A. and Kendall, W.S. (2002)
  Simulation of cluster point processes without edge effects.
  \emph{Advances in Applied Probability} \bold{34}, 267--280.
}
\seealso{
  \code{\link[spatstat.random]{rNeymanScott}},
  \code{\link[spatstat.random]{rMatClust}},
  \code{\link[spatstat.random]{rThomas}},
  \code{\link[spatstat.random]{rCauchy}},
  \code{\link[spatstat.random]{rVarGamma}}
}
\examples{
  Y <- rclusterBKBC("Thomas", 10,5,0.2)
  Y
  Z <- rclusterBKBC("VarGamma", 10,5,0.2,
          nu=-1/4,
          internal="naive", external="super",
          verbose=FALSE)
}
\keyword{datagen}
\keyword{spatial}
