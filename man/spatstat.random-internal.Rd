\name{spatstat.random-internal} 
\title{Internal spatstat.random functions}
\alias{spatstat.random-internal} %DoNotExport

\alias{change.default.expand}          
\alias{datagen.runifpointOnLines}
\alias{datagen.runifpoisppOnLines}
\alias{datagen.rpoisppOnLines}
\alias{default.clipwindow}
\alias{detect.par.format}
\alias{expandwinPerfect}
\alias{fakeNeyScot}
\alias{getRandomFieldsModelGen}
\alias{handle.rshift.args}
\alias{HermiteCoefs}
\alias{is.cadlag}
\alias{is.expandable}
\alias{is.expandable.rmhmodel}
\alias{kraever}
\alias{kraeverRandomFields}
\alias{MultiPair.checkmatrix}
\alias{optimalinflation}
\alias{print.rmhcontrol}
\alias{print.rmhexpand}
\alias{print.rmhmodel}
\alias{print.rmhstart}
\alias{print.rmhInfoList}
\alias{print.summary.rmhexpand}
\alias{RandomFieldsSafe}
\alias{reheat}
\alias{resolve.vargamma.shape}
\alias{retrieve.param}
\alias{rMaternInhibition}
\alias{rMatClustHom}
\alias{rThomasHom}
\alias{rCauchyHom}
\alias{rGRFcircembed}
\alias{rGRFgauss}
\alias{rGRFexpo}
\alias{rGRFstable}
\alias{rGRFgencauchy}
\alias{rGRFmatern}
\alias{rPoissonClusterEngine}
\alias{RmhExpandRule}
\alias{rmhsnoop}
\alias{rmhResolveControl}
\alias{rmhResolveExpansion}
\alias{rmhResolveTypes}
\alias{rmhSnoopEnv}
\alias{rmhcontrol.rmhcontrol}
\alias{rmhcontrol.list}
\alias{rmhEngine}
\alias{rmhmodel.rmhmodel}
\alias{rmhstart.rmhstart}
\alias{rmhstart.list}
\alias{rmpoint.I.allim}
\alias{rpoint.multi}
\alias{rthinEngine}
\alias{runifpoispp}          
\alias{runifpoisppOnLines}          
\alias{spatstatClusterModelInfo}
\alias{spatstatClusterSimInfo}
\alias{spatstatClusterSimModelMatch}
\alias{spatstatRmhInfo}
\alias{summarise.trend}
\alias{summary.rmhexpand}
\alias{thinjump}
\alias{thinParents}
\alias{update.rmhstart}
\alias{validate.kappa.mu}



%%%%%%%
\description{
  Internal spatstat.random functions.
}
\usage{



change.default.expand(x, newdefault)
datagen.runifpointOnLines(n, L)
datagen.runifpoisppOnLines(lambda, L)
datagen.rpoisppOnLines(lambda, L, lmax, \dots, check)
default.clipwindow(object, epsilon)
expandwinPerfect(W, expand, amount)
detect.par.format(par, native, generic)
fakeNeyScot(Y, lambda, win, saveLambda, saveparents)
getRandomFieldsModelGen(model)
handle.rshift.args(W, \dots, radius, width, height,
                   edge, clip, edgedefault)
HermiteCoefs(order)
is.cadlag(s)
is.expandable(x)
\method{is.expandable}{rmhmodel}(x) 
kraever(package, fatal)
kraeverRandomFields()
MultiPair.checkmatrix(mat, n, matname, naok, zerook, asymmok)
optimalinflation(clusters, mod, rD)
\method{print}{rmhcontrol}(x, \dots)
\method{print}{rmhexpand}(x, \dots, prefix=TRUE)
\method{print}{rmhmodel}(x, \dots)
\method{print}{rmhstart}(x, \dots)
\method{print}{rmhInfoList}(x, \dots)
\method{print}{summary.rmhexpand}(x, \dots)
RandomFieldsSafe()
reheat(model, invtemp)
resolve.vargamma.shape(\dots, nu.ker, nu.pcf, nu, allow.nu, allow.default)
retrieve.param(desired, aliases, \dots, par)
rMaternInhibition(type, kappa, r, win, stationary, \dots, nsim, drop)
rMatClustHom(kappa, mu, R, W, \dots, nsim, drop, inflate, saveparents)
rThomasHom(kappa, mu, sigma, W, \dots,
           nsim, drop, inflate, saveparents, maxinflate)
rCauchyHom(kappa, mu, scale, W, \dots,
           nsim, drop, inflate, saveparents, maxinflate)
rGRFcircembed(W, mu, var, corrfun, \dots, nsim, drop)
rGRFgauss(W, mu, var, scale, \dots, nsim, drop)
rGRFexpo(W, mu, var, scale, \dots, nsim, drop)
rGRFstable(W, mu, var, scale, alpha, \dots, nsim, drop)
rGRFgencauchy(W, mu, var, scale, alpha, beta, \dots, nsim, drop)
rGRFmatern(W, mu, var, scale, nu, \dots, nsim, drop)
rPoissonClusterEngine(kappa, expand, rcluster, win, \dots,
                      nsim, drop, saveparents, kappamax, lmax, rmax) 
RmhExpandRule(nama)
\method{rmhcontrol}{rmhcontrol}(\dots) 
\method{rmhcontrol}{list}(\dots) 
rmhEngine(InfoList, \dots, verbose, kitchensink, preponly, snoop,
                           overrideXstart, overrideclip) 
rmhResolveControl(control, model)
rmhResolveExpansion(win, control, imagelist, itype)
rmhResolveTypes(model, start, control)
rmhsnoop(\dots, Wsim, Wclip, R, xcoords, ycoords, mlevels, mcodes, irep, itype, 
     proptype, proplocn, propmark, propindx, numerator, denominator,
     panel.only)
rmhSnoopEnv(Xinit, Wclip, R)
\method{rmhmodel}{rmhmodel}(model, \dots) 
\method{rmhstart}{rmhstart}(start, \dots) 
\method{rmhstart}{list}(start, \dots) 
rmpoint.I.allim(n, f, types, nsim)
rpoint.multi(n, f, fmax, marks, win, giveup, verbose, warn, nsim, drop)
rthinEngine(X, P, \dots, nsim, drop, Pmax, na.zero, what, fatal, warn)
runifpoispp(lambda, win, \dots, nsim, drop)
runifpoisppOnLines(lambda, L, nsim, drop)
spatstatClusterModelInfo(name, onlyPCP)
spatstatClusterSimInfo(name)
spatstatClusterSimModelMatch(name, verbose)
spatstatRmhInfo(cifname)
summarise.trend(trend, w, a)
\method{summary}{rmhexpand}(object, \dots)
thinjump(n, p)
thinParents(X, P, Pmax)
\method{update}{rmhstart}(object, \dots)
validate.kappa.mu(kappa, mu, kappamax, mumax,
                  win, expand, \dots, context)


}
\details{
  These internal \pkg{spatstat.random} functions should not be called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat.random} to the next.
}
\value{
The return values of these functions are not documented,
and may change without warning.
}
\keyword{internal}
