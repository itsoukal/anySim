#' @title Estimation of the auxiliary n-AR(1) model parameters
#'
#' @description Estimation of parameters of the the sum of n univariate AR(1) models to simulate the auxiliary Gaussian process.
#'
#' @param ACF A vector with the target autocorrelation structure (including lag-0, i.e., 1).
#' @param Ar1Num A scalar (>=2) indicating the number (n) of AR(1) models.
#' @param dist A string indicating the quantile function of the target marginal distribution (i.e., the ICDF).
#' @param params A named list with the parameters of the target distribution.
#' @param NatafIntMethod A string ("GH", "Int", or "MC"), indicating the intergation method, to resolve the Nataf integral.
#' @param NoEval A scalar indicating (default: 9) the number of evaluation points for the integration methods.
#' @param polydeg A scalar indicating the order of the fitted polynomial. If polydeg=0, then another curve is fitted.
#' @param ... Additional named arguments for the selected "NatafIntMethod" method.
#'
#' @note Avoid the use of the "GH" method (i.e., NatafIntMethod='GH'), when the marginal(s) are discrete.
#'
#' @return A list with the parameters of the auxiliary Gaussian n-AR(1) model.
#' @export
#'
#' @examples
#' ## Parameter estimation for a process with zero-inflated (i.e., mixed) marginal distribution,
#' ## with p0=0.8, and a Gamma distribution
#' ## for the continuous part with shape=0.5 and scale=1.
#' ## In this case, the Autocorrelation strucure is a CAS ACS with b=2 and k=0.5.
#'\dontrun{
#' ACF=acsCAS(param = c(2, 0.5), lag = 500, var = 1)
#' dist='qmixed'
#' params=list(Distr=qgamma, p0=0.8, shape=0.5, scale=1)
#' nAR1param=EstnAR1(ACF = ACF, Ar1Num = 3, dist = dist, params = params,
#'                   NatafIntMethod = 'GH', NoEval = 9, polydeg = 8)
#'}
EstnAR1<-function(ACF, Ar1Num, dist='qgamma', params, NatafIntMethod = 'GH', NoEval = 9, polydeg = 8, ...){
  Nataf=NatafInvD(targetrho = ACF[-1], fx = dist, fy = dist,
                  paramlistfx = params, paramlistfy = params,
                  NatafIntMethod = NatafIntMethod,
                  NoEval = NoEval, polydeg = polydeg, ...)

  GAFv=(c(1,Nataf$rzEq))
  Natafdf=Nataf$dfnataf

  lbr=rep(0.0001,Ar1Num)
  ubr=rep(0.9999,Ar1Num)
  lbs=rep(0.0001,Ar1Num)
  ubs=rep(0.9999,Ar1Num)

  lb=c(lbr,lbs)
  ub=c(ubr,ubs)

  x0=runif(Ar1Num*2) # initial values

  #TO-DO: drop the depencence in Rsolnp

  Res=solnp(pars = x0, fun = nAR1obj, eqfun = nAR1con,
            eqB=1, LB=lb, UB=ub,
            ACF=GAFv, Ar1Num=Ar1Num)

  optPars=Res$pars
  optValue=nAR1obj(optPars,GAFv,Ar1Num)
  return(list('oPars'=optPars, 'oValue'=optValue, 'GAF'=GAFv, 'dist'=dist, 'params'=params, 'Natafdf'=Natafdf, 'Ar1Num'=Ar1Num))
}

nAR1obj<-function(par,ACF,Ar1Num){
  len=length(ACF[-1])
  lag=seq(1:len)

  r=par[1:Ar1Num]
  s=par[(Ar1Num+1):length(par)]

  X=matrix(NA,len,Ar1Num)

  for (i in 1:Ar1Num){
    for (j in 1:len) {
      X[j,i]=s[i]*r[i]^j

    }
  }
  XX=apply(X,1,sum)
  sse=sum((ACF[-1]-XX)^2)
  return(sse)
}

nAR1con<-function(par,ACF,Ar1Num) {
  c=sum(par[(Ar1Num+1):length(par)])
  return(c)
}
