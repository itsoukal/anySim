#' @title Estimation of parameters of the the auxiliary Gaussian AR(p) model.
#'
#' @description Estimation of parameters of AR(p) model to simulate the auxiliary Gaussian process.
#'
#' @param ACF A vector with the target autocorrelation structure (including lag-0 coefficient that is equal to 1).
#' @param maxlag A scalar incating the order of the AR(p) model. If maxlag=0, then the order of the model is p=(length(ACF)-1)
#' @param dist A string indicating the quantile function of the target marginal distribution (i.e., the ICDF).
#' @param params A named list with the parameters of the target distribution.
#' @param NatafIntMethod A string ("GH", "Int", or "MC") indicating the intergation method to resolve the Nataf integral.
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
#' ## Simulation of univariate stationary process with Gamma marginal distribution 
#' ## and autocorrelation structure given by the product of a CAS and a periodic ACS.
#'\dontrun{
#' set.seed(12)
#' 
#' # Define the target autocorrelation structure.
#' acsS=csCAS(param=c(3,0.6),lag=1000) # Stationary CAS with b=3 and k=0.6.
#' acsP=csPeriodic(param=c(12,1.5),lag=1000) # Periodic ACS with p=12 and l=1.5.
#' ACS=csP*csS # The target ACS as product of the two previous ones.
#' 
#' # Define the target distribution function (ICDF).
#' FX='qgamma' # the Gamma distribution
#' 
#' # Define the parameters of the target distribution.
#' pFX=list(shape=5,scale=1)
#' 
#' # Estimate the parameters of the auxiliary Gaussian AR(p) model.
#' ARTApar=EstARTAp(ACF=ACS,maxlag=0,dist=FX,params=pFX,NatafIntMethod='GH')
#' 
#' # Generate a synthetic series of 10000 length. 
#' SynthARTAcont=SimARTAp(ARTApar=ARTApar,steps=10^5)
#'}
#'
#' ## Simulation of univariate stationary process with discrete marginal distribution
#' ## (Beta-Binomial) and autocorrelation structure given by CAS.
#'\dontrun{
#' set.seed(16)
#' 
#' # Define the target autocorrelation structure.
#' ACS=acsCAS(param=c(1.5,0.3),lag=1000) # CAS with b=1.5 and k=0.3.
#' 
#' # Define the target distribution function (ICDF).
#' require(TailRank)
#' FX='qbb' # the Beta-Binomial distribution.
#' 
#' # Define the parameters of the target distribution.
#' pFX=list(N=10,u=3,v=10)
#' 
#' # Estimate the parameters of the auxiliary Gaussian AR(p) model.
#' ARTApar=EstARTAp(ACF=ACS,maxlag=0,dist=FX,params=pFX,NatafIntMethod="MC")
#' 
#' # Generate a synthetic series of 10000 length. 
#' SynthARTAdiscr=SimARTAp(ARTApar=ARTApar,steps=10^5)
#'}
#'
#' ## Simulation of univariate stationary process with zero-inflated marginal distribution 
#' ## (Gen. Gamma for the continuous part) and autocorrelation structure given by CAS.
#'\dontrun{
#' set.seed(18)
#' 
#' # Define the target autocorrelation structure.
#' ACS=acsCAS(param=c(0.91,1.09),lag=1000) # CAS with b=0.91 and k=1.09.
#' 
#' # Define the target distribution function (ICDF).
#' FX='qmixed' # Define that distribution is of zero-inflated type.
#' 
#' # Define the distribution for the continuous part of the process.
#' # Here, a re-parameterized version of Gen. Gamma distribution is used.
#' qgengamma=function(p,scale,shape1,shape2){
#'   require(VGAM)
#'   X=qgengamma.stacy(p=p,scale=scale,k=(shape1/shape2),d=shape2)
#'   return(X)
#' } 
#' 
#' # Define the parameters of the zero-inflated distribution function.
#' pFX=list(Distr=qgengamma,p0=0.8,scale=0.25,shape1=1.16,shape2=0.54) 
#' 
#' # Estimate the parameters of the auxiliary Gaussian AR(p) model.
#' ARTApar=EstARTAp(ACF=ACS,dist=FX,params=pFX,NatafIntMethod="GH",NoEval=9,polydeg=0)
#' 
#' # Generate a synthetic series of 10000 length. 
#' SynthARTAzi=SimARTAp(ARTApar=ARTApar,steps=10^5)
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
