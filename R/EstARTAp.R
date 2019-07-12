#' @title Estimation of the auxiliary AR(p) model parameters
#'
#' @description Estimation of parameters of the auxiliary AR(p) model to simulate the auxiliary Gaussian process.
#'
#' @param ACF A vector with the target autocorrelation structure (including lag-0, i.e., 1).
#' @param maxlag A scalar incating the order of the AR(p) model. If maxlag=0, then the order of the model is p=(length(ACF)-1)
#' @param dist A string indicating the quantile function of the target marginal distribution (i.e., the ICDF).
#' @param params A named list with the parameters of the target distribution.
#' @param NatafIntMethod A string ("GH", "Int", or "MC"), indicating the intergation method, to resolve the Nataf integral.
#' @param NoEval A scalar indicating (default: 9) the number of evaluation points for the integration methods.
#' @param polydeg A scalar indicating the order of the fitted polynomial. If polydeg=0, then another curve is fitted.
#' @param ... Additional named arguments for the selected "NatafIntMethod" method.
#'
#' @note Avoid the use of the "GH" method (i.e., NatafIntMethod='GH'), when the marginal(s) are discrete.
#'
#' @return A list with the parameters of the auxiliary Gaussian AR(p) model.
#' @export
#'
#' @examples
#' ## Parameter estimation for a process with zero-inflated (i.e., mixed) marginal distribution,
#' ## with p0=0.9, and a Gamma distribution
#' ## for the continuous part with shape=0.1 and scale=1.
#' ## In this case, the Autocorrelation strucure is a simple AR(1) with rho=0.8.
#'\dontrun{
#' ACF=as.vector(ARMAacf(ar = 0.8, lag.max = 100))
#' fx='qmixed'
#' pfx=list(Distr=qgamma, p0=0.9, shape=0.1, scale=1)
#'
#' ARTApar=EstARTAp(ACF=ACF, maxlag=0, dist=fx, params=pfx,
#' NatafIntMethod ='GH', NoEval=9, polydeg=0)
#'}
EstARTAp=function(ACF, maxlag=0, dist, params, NatafIntMethod='GH', NoEval=9, polydeg=8, ...){

  N=length(ACF)
  p=N-1
  if (maxlag==0){ maxlag=p }

  if (NatafIntMethod=='GH') {
    Nataf=NatafInvD(targetrho = ACF[-1], fx = dist, fy = dist, paramlistfx = params, paramlistfy = params,
                    NatafIntMethod = 'GH', polydeg=polydeg, NoEval = NoEval, ...)
  } else if (NatafIntMethod=='Int') {
    Nataf=NatafInvD(targetrho = ACF[-1], fx = dist, fy = dist, paramlistfx = params, paramlistfy = params,
                    NatafIntMethod = 'Int', polydeg=polydeg, NoEval = NoEval)
  } else if (NatafIntMethod=='MC') {
    Nataf=NatafInvD(targetrho = ACF[-1], fx = dist, fy = dist, paramlistfx = params, paramlistfy = params,
                    NatafIntMethod = 'MC', polydeg=polydeg, NoEval = NoEval, ...)
  } else {
    print('Select an appropriate method to resolve Nataf integral. That is: "GH, "MC" or "Int"')
  }

  Natafdf=Nataf$dfnataf
  ACFn=Nataf$rzEq
  ACFn=ifelse(is.na(ACFn),0,ACFn)
  ACFn=c(1,ACFn)
  # Solve the Yule-Walker system
  R=toeplitz(ACFn)[1:maxlag,1:maxlag]
  phi=solve(R,ACFn[-1][1:maxlag])
  sigma=sqrt(1-(ACFn[-1]%*%phi))
  # Store the parameters into a list.
  ARpar=list('ACFn'=ACFn,'phi'=phi,'simga'=sigma,'dist'=dist,'params'=params,'Natafdf'=Natafdf)
  return(ARpar)
}
