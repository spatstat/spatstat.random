##  spatstat.random/R/First.R

.onLoad <- function(...) reset.spatstat.options()

.onAttach <- function(libname, pkgname) {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat.random"),
                 fields="Version")
  vs <- as.character(vs)
  putSpatstatVariable("SpatstatRandomVersion", vs)
  packageStartupMessage(paste("spatstat.random", vs))
  return(invisible(NULL))
}

  
