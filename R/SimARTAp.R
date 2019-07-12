#' @title Simulation of the target stationary process using the ARTA(p) model.
#'
#' @description Simulation of the target stationary process using an ARTA(p) model to simulate the auxiliary Gaussian process and establish the target correlation structure.
#'
#' @param ARTApar A list containing the parameters of the model. The list is constructed by the function "EstARTAp".
#' @param burn A scalar specifying the length of burn-out sample.
#' @param steps A scalar specifying the length of the time series to be generated.
#' @param stand A boolean (T or F) indicating whether to standardize (or not) the auxiliary Gaussian time series prior to their mapping to the actual domain. The default value is FALSE.
#'
#' @return A list of the 3 generated time series (in vector format):
#' X: The final time series at the actual domain with the target marginal distribution and correlation structure;
#' Z: The auxiliary Gaussian time series at the Gaussian domain and;
#' U: The auxiliary uniform time series at the Copula domain (i.e., in [0,1]).
#'
#' @export
#'
#' @examples
#'
#' ## Simulation of a process with a zero-inflated (i.e., mixed)
#' ## marginal distribution, with p0=0.9, and a Gamma distribution
#' ## for the continuous part with shape=0.1 and scale=1.
#' ## In this case, the target autocorrelation strucure is from
#' ## the CAS ACS with b=1 and k=0.6.
#' \dontrun{
#' ACF=acsCAS(param=c(1, 0.6), lag=100, var=1)
#' fx='qmixed'
#' pfx=list(Distr=qgamma, p0=0.9, shape=0.1, scale=1)
#'
#' ARTApar=EstARTAp(ACF=ACF, maxlag=0, dist=fx, params=pfx,
#' NatafIntMethod ='GH', NoEval=9, polydeg=8)
#'
#' Sim=SimARTAp(ARTApar = ARTApar, burn = 1000, steps = 10^5, stand = 0)
#' acf(Sim$X)
#' lines(0:(length(ACF)-1), ACF)
#' plot(Sim$X[1:1000], type='l', col='red')
#'
#' ## Simulation of a process with a bernoulli marginal distribution,
#' ## with size=1, and prob=0.2.
#' ## In this case, the target autocorrelation strucure
#' ## is from an fGn process (i.e., Hurst) with H=0.7.
#'
#' ACF=acsHurst(H=0.7, lag=500, var=1)
#' fx='qbinom'
#' pfx=list(size=1, prob=0.2)
#'
#' ARTApar=EstARTAp(ACF=ACF, maxlag=0, dist=fx, params=pfx,
#' NatafIntMethod ='Int', NoEval=9, polydeg=0)
#'
#' Sim=SimARTAp(ARTApar = ARTApar, burn = 1000, steps = 10^5, stand = 0)
#' acf(Sim$X)
#' lines(0:(length(ACF)-1), ACF)
#' plot(Sim$X[1:1000], type='l', col='red')
#'
#' ## Simulation of a process with a Beta marginal distribution,
#' ## with shape=2, and shape2=10.
#' ## In this case, the target autocorrelation strucure is from
#' ## an fGn process (i.e., Hurst) with H=0.7.
#'
#' ACF=acsHurst(H=0.7, lag=500, var=1)
#' fx='qbeta'
#' pfx=list(shape1=2, shape2=10)
#'
#' ARTApar=EstARTAp(ACF=ACF, maxlag=0, dist=fx, params=pfx,
#' NatafIntMethod ='Int', NoEval=9, polydeg=0)
#'
#' Sim=SimARTAp(ARTApar = ARTApar, burn = 1000, steps = 10^5, stand = 0)
#' acf(Sim$X)
#' lines(0:(length(ACF)-1), ACF)
#' plot(Sim$X[1:1000], type='l', col='red')
#'}
SimARTAp=function(ARTApar, burn=1000, steps=10000, stand=F){
  phi=ARTApar$phi
  sigma=ARTApar$simga
  dist=ARTApar$dist
  params=ARTApar$params

  lag=length(phi)
  Start=rnorm(lag,mean = 0,sd = 1)
  Z=mat.or.vec((burn+steps),1)
  Z[1:lag]=Start

  totalSteps=burn+steps
  for (i in (lag+1):totalSteps){
    Z[i]=Z[i-c(1:lag)]%*%phi+rnorm(1,mean = 0,sd=sigma)
  }
  Z=Z[(burn+1):totalSteps]

  if (stand==T){
    U=pnorm(Z,mean = mean(Z), sd=sd(Z))
    Z=qnorm(U)
  } else {
    U=pnorm(Z)
  }

  funcall <- as.call(c(as.name(dist), list(U), params))
  X=eval(funcall)

  return(list('X'=X,'Z'=Z,'U'=U))
}
