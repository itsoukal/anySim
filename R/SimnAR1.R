#' @title Simulation of the target stationary process using the n-ARTA(1) model.
#'
#' @description Simulation of the target stationary process using an n-ARTA(1) model to simulate an auxiliary Gaussian process and establish the target correlation structure.
#'
#' @param nAR1param A list containing the parameters of the model. The list is constructed by the function "estnAR1".
#' @param steps A scalar specifying the length of the time series to be generated.
#'
#' @return A list with 3 generated time series (in vector format):
#' X: The final time series at the actual domain with the target marginal distribution and correlation structure;
#' Z: The auxiliary Gaussian time series at the Gaussian domain and;
#' U: The auxiliary uniform time series at the Copula domain (i.e., in [0,1]).
#'
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
#'
#' Sim=SimnAR1(nAR1param = nAR1param, steps = 10^6)
#'
#' acf(Sim$X)
#' lines(0:500,ACF,col='red')
#' plot.ecdf(Sim$X)
#'}
#'
SimnAR1<-function(nAR1param, steps) {
  par=nAR1param$oPars
  dist=nAR1param$dist
  params=nAR1param$params
  Ar1Num=nAR1param$Ar1Num

  r=par[1:Ar1Num]
  s=par[(Ar1Num+1):length(par)]

  r2=r^2
  s2=sqrt((1-r2)*s)


  X=matrix(NA,steps,Ar1Num)
  X[1,]=rnorm(Ar1Num)

  for (i in 1:Ar1Num) {
    for (j in 2:steps){
      X[j,i]=r[i]*X[j-1,i]+rnorm(1,0,sd=s2[i])
    }
  }
  XX=apply(X,1,sum)
  XX=matrix(XX)
  U=pnorm(XX)
  X=eval(as.call(c(as.name(dist), list(U), params)))
  return(list('X'=X, 'Z'=XX, 'U'=U))
}
# ACF=acsCAS(param = c(2, 0.5), lag = 500, var = 1)
# dist='qmixed'
# params=list(Distr=qgamma, p0=0.8, shape=0.5, scale=1)
# nAR1param=EstnAR1(ACF = ACF, Ar1Num = 3, dist = dist, params = params,
#                   NatafIntMethod = 'GH', NoEval = 9, polydeg = 8)
#
# Sim=SimnAR1(nAR1param = nAR1param, steps = 10^6)
#
# acf(Sim$X)
# lines(0:500,ACF,col='red')
