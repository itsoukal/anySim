#' @title Zero-Inflated distribution model (i.e., mixed)
#'
#' @description Density, distribution function, quantile function and random generation
#'  for the zero-inflated (i.e., mixed) distribution model. This model is composed by two parts,
#'  1) the discrete part, which regards an atom at zero, and
#'  2) the continuous part, which regards a continuous distribution model (J-shaped and with left support >0).
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param Distr The name (as a function) of the continuous distribution model.
#' @param p0 Probability of zero values (i.e., zero-inflation).
#' @param ... Additional named arguments containing the continuous distribution parameters
#'
#' @name mixed
#' @aliases mixed
#' @aliases dmixed
#'
#' @return
#' dmixed gives the density of the zero-inflated (i.e., mixed) distribution model.
#' pmixed gives the cdf of the zero-inflated (i.e., mixed) distribution model.
#' qmixed gives the quantile (ICDF) of the zero-inflated (i.e., mixed) distribution model.
#' rmixed gives random variates from the zero-inflated (i.e., mixed) distribution model.
#' @export
#'
#' @examples
#' ## Plot the CDF of a Gamma distribution.
#' p=seq(0,1,0.01)
#' x=qmixed(p, Distr = qgamma, p0=0.7, shape=0.5, scale=1)
#' plot(x, p)
#'
#' ## Generate 100000 random variables with p0=0.7 and
#' ## Gamma distribution for the continoous part.
#' X=rmixed(1000, qgamma, p0=0.7, shape=0.5, scale=1)
#' hist(X)
#'
#' ## Generate 100000 random variables with p0=0.7 and
#' ## Burr type XII distribution for the continuous part.
#' ## The actuar package is required, since it contains
#' ## the d,p,q,r functions of the Burr type XII distribution.
#' require(actuar)
#' X=rmixed(1000, qburr, p0=0.7, shape1=5, shape2=1, scale=1)
#' hist(X)
#' plot(sort(X[X>0]), 1-ppoints(X[X>0], a=0), log='xy',
#' xlab = 'x', ylab = 'P[X>x]', main='Prob of exceedance plot')
#'
dmixed=function(x,Distr,p0,...) {
  y=rep(NA,length(x))

  for (i in 1:length(x)){
  if (x[i]==0){
    y[i]=p0
  }else {
    y[i]=(1-p0)*Distr(x[i],...)
  }
  }

  return(y)
}

#' @rdname mixed
#' @export

pmixed=function(q, Distr, p0,...) {
  y=rep(NA,length(q))

  for (i in 1:length(q)){
    y[i]=(1-p0)*Distr(q[i],...)+p0

  }

  return(y)
}

#' @rdname mixed
#' @export

qmixed=function(p, Distr, p0,...) {
  y=rep(NA,length(p))

  for (i in 1:length(p)){
    if (p[i]<=p0){
      y[i]=0
    }else {
      u=(p[i]-p0)/(1-p0)
      y[i]=Distr(u,...)
    }
  }

  return(y)
}

#' @rdname mixed
#' @export

rmixed=function(n, Distr, p0, ...) {
  y=qmixed(p = runif(n), Distr = Distr, p0 = p0, ...)
  return(y)
}
