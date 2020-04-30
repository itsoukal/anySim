#' @title Solve the Nataf integral with an alternative integration method
#'
#' @description Estimation of the resulting (i.e., in the actual domain) correlation coefficients,
#' given the equivalent correlation coefficients (i.e., in the Gaussian domain).
#'
#' @param rho A scalar or vector of correlation coefficients (i.e., in seq(from=0, to=1, by=0.1)).
#' @param fx A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param fy A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param paramlistfx A named list with the parameters of the distribution.
#' @param paramlistfy A named list with parameters of the distribution.
#'
#' @return A vector of correlation coefficients in the actual domain.
#' @export
#'
#' @examples
#' ## The case of two identrical Gamma distributions, with shape=1 and scale=1.
#'\dontrun{
#' fx=fy='qgamma'
#' pfx=pfy=list(shape=1, scale=1)
#' rhoz=seq(from=0, to=1 , by=0.2)
#' rhox=NatafInt(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' plot(rhoz,rhox); abline(0,1)
#'
#' ## The case identrical Bernoulli distributions, with size=1 and prob=0.3.
#' fx=fy='qbinom'
#' pfx=pfy=list(size=1, prob=0.3)
#' rhoz=seq(from=0, to=1 , by=0.2)
#' rhox=NatafInt(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' plot(rhoz,rhox); abline(0,1)
#'
#' ## The case of two identrical zero-inflated (i.e., mixed) distributions,
#' ## with p0=0.7 a Gamma distribution
#' ## for the continuous part with shape=1 and scale=1.
#'
#' fx=fy='qzi'
#' pfx=pfy=list(Distr=qgamma, p0=0.7, shape=1, scale=1)
#' rhoz=seq(from=0, to=1 , by=0.2)
#' rhox=NatafInt(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' plot(rhoz,rhox);abline(0,1)
#' # compare with the Gauss-Hermite integration method
#' rhox=NatafGH(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy, nodes = 21)
#' points(rhoz,rhox,col='red',pch=19)
#'}
NatafInt=function(rho, fx, fy, paramlistfx, paramlistfy){

  if (fx !='qzi') {
    funfx=match.fun(fx)
    stats1=do.call(DistrStats2, args = c(qDistr=funfx, paramlistfx))

  } else {
    funfx=match.fun(fx)
    p0=paramlistfx$p0
    Distr=paramlistfx$Distr
    pfx=paramlistfx
    pfx$Distr=NULL
    pfx$p0=NULL
    suppressWarnings( stats1<-DistrStats2(qDistr = function(x) funfx(x, Distr = Distr, p0 = p0,  unlist(pfx))) )
  }
  if (fy !='qzi') {
    funfy=match.fun(fy)
    stats2=do.call(DistrStats2, args = c(qDistr=funfy, paramlistfy))
  } else {
    funfy=match.fun(fy)
    p0=paramlistfy$p0
    Distr=paramlistfy$Distr
    pfy=paramlistfy
    pfy$Distr=NULL
    pfy$p0=NULL
    suppressWarnings( stats2<-DistrStats2(qDistr = function(x) funfy(x, Distr = Distr, p0 = p0,  unlist(pfy))) )
  }

  mu1=stats1[1]
  sigma1=stats1[2]^0.5
  mu2=stats2[1]
  sigma2=stats2[2]^0.5
  Q1=-(mu1*mu2)/(sigma1*sigma2)
  Q2=1/(2*pi*sigma1*sigma2)

  Rx=rho
  lim=7
  substr(fx,start = 1,stop = 1)='q'
  substr(fy,start = 1,stop = 1)='q'

  for (i in 1:length(rho)){

    Intres=cubature::adaptIntegrate(f = Nataf2dIntegral,lowerLimit = rep(-lim,2),upperLimit = rep(lim,2),
                                    maxEval = 10000,
                                    Rho=rho[i], fx=fx, fy=fy,
                                    paramlistfx=paramlistfx, paramlistfy=paramlistfy)


    while (is.infinite(Intres$integral) || is.nan(Intres$integral) ){
      lim=lim-1
      Intres=cubature::adaptIntegrate(f = Nataf2dIntegral,lowerLimit = rep(-lim,2),upperLimit = rep(lim,2),
                                      maxEval = 10000,
                                      Rho=rho[i], fx=fx, fy=fy,
                                      paramlistfx=paramlistfx, paramlistfy=paramlistfy)
    }
    Rx[i]=Q1+Q2*Intres$integral
    lim=7
  }

  return(Rx)
}

Nataf2dIntegral=function(U,Rho,fx,fy,paramlistfx,paramlistfy){
  Int=(eval((as.call(c(as.name(fx), pnorm(U[1]), paramlistfx)))))*
    (eval((as.call(c(as.name(fy), pnorm((Rho*U[1]+sqrt(1-Rho^2)*U[2])), paramlistfy)))))*
    (exp(-((U[1]^2)+(U[2]^2))/2))
  return(Int)
}

