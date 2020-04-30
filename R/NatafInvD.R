#' @title Direct estimation of equivalent correlation coefficients.
#'
#' @description Direct estimation of equivalent correlation coefficients (i.e., in the Gaussian domain).
#'
#' @param targetrho A scalar or vector of target correlation coefficients.
#' @param fx A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param fy A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param paramlistfx A named list with the parameters of the distribution.
#' @param paramlistfy A named list with parameters of the distribution.
#' @param NoEval A scalar indicating (default: 9) the number of evaluation points for the integration methods.
#' @param NatafIntMethod A string ("GH", "Int", or "MC"), indicating the integration method, to resolve the Nataf integral.
#' @param polydeg A scalar indicating the order of the fitted polynomial. If polydeg=0, then an alternative parametric curve is fitted (see, Papalexiou, 2018).
#' @param ... Additional named arguments for the selected "NatafIntMethod" method.
#'
#' @note Avoid the use of the "GH" method (i.e., NatafIntMethod='GH'), when the marginal(s) are discrete.
#' @return A named list with two elements:
#' dfnataf: A dataframe that contains the pairs of Gaussian and resulting correlation coefficients, upon which the curve (polynomial or other) was fitted.
#' rzEq: A vector with the equivalent correlation coefficients, that result into the target ones (i.e., targetrho).
#'
#' @export
#'
#' @examples
#' ## The case of two identrical zero-inflated (i.e., mixed) distributions,
#' ## with p0=0.9 a Gamma distribution
#' ## for the continuous part with shape=0.1 and scale=1.
#'\dontrun{
#' fx=fy='qzi'
#' pfx=pfy=list(Distr=qgamma, p0=0.9, shape=0.1, scale=1)
#' rhoz=seq(from=0, to=1 , length.out = 21)
#' rhox=NatafInt(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' plot(rhox,rhoz, col='red', pch=19); abline(0,1)
#' rhotarget=seq(from=0.0001, to=0.9999 , length.out = 210)
#' req=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#'               NatafIntMethod = 'GH', polydeg=8, NoEval = 9)
#' points(rhotarget, req, col='blue', pch=17);
#'
#' req=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#'               NatafIntMethod = 'GH', polydeg=8, NoEval = 9)
#' points(rhotarget, req, col='green')
#'}
#' ## The case with identical Bernoulli distributions, with size=1 and prob=0.2.
#'
#'\dontrun{
#' fx=fy='qbinom'
#' pfx=pfy=list(size=1, prob=0.2)
#' rhoz=seq(from=0, to=1 , length.out = 21)
#' rhox=NatafInt(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' plot(rhox,rhoz, col='red', pch=19); abline(0,1)
#' rhotarget=seq(from=0.0001, to=0.9999 , length.out = 210)
#' req=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#'               NatafIntMethod = 'Int', polydeg=8, NoEval = 9)$rzEq
#' points(rhotarget, req, col='blue', pch=17);
#'
#' req2=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#'               NatafIntMethod = 'Int', polydeg=8, NoEval = 9)$rzEq
#' points(rhotarget, req2, col='green')
#'}
#'
NatafInvD=function(targetrho, fx, fy, paramlistfx, paramlistfy, NatafIntMethod='GH', NoEval=19, polydeg=8, ...){

  rmax=ifelse(max(targetrho)>0,1,0)
  rmin=ifelse(min(targetrho)<0,-1,0)

  rz=seq(rmin, rmax, length.out = NoEval)
  Index0=which(rz==0)

  if (NatafIntMethod=='GH') {
    rx=NatafGH(rho = rz, fx = fx, fy = fy, paramlistfx = paramlistfx, paramlistfy = paramlistfy, ...)
  } else if (NatafIntMethod=='MC') {
    rx=NatafMC(rho = rz, fx = fx, fy = fy, paramlistfx = paramlistfx, paramlistfy = paramlistfy, ...)
  } else if (NatafIntMethod=='Int') {
    rx=NatafInt(rho = rz, fx = fx, fy = fy, paramlistfx = paramlistfx, paramlistfy = paramlistfy)
  } else {
    print('Select an appropriate method to resolve Nataf integral. That is: "GH, "MC" or "Int"')
  }
  rx[Index0]=0
  dfnataf=data.frame(rz=rz, rx=rx)


  if (polydeg!=0){
    lb=min(rz);
    ub=max(rz)
    x=seq(lb, ub, 0.0002)
    fitZ2Y=pracma::polyfit(x = rz, y = rx, n = polydeg)
    y=pracma::polyval(fitZ2Y,x)
    X=(cbind(y,x))
    rzEq=approx(X[,1], X[,2], targetrho)$y
  }else {
    rxmax=rx[length(rx)]
    rxmax=ifelse(rxmax>1,1,rxmax)

    DE=DEoptim::DEoptim(fn = TwoParfunObj, lower = c(0.0001,0.0001),upper = c(10000, 10000),
                        rx=rx, rz=rz, rmax=rxmax,
                        control = DEoptim.control(trace = 0,NP=30,itermax =1000))
    DE$optim$bestmem
    b=DE$optim$bestmem[1]
    c=DE$optim$bestmem[2]
    targetrho=ifelse(targetrho>rxmax,NA,targetrho)
    rzEq=TwoParfun(b, c, rx=targetrho, rmax = rxmax)
    #
    # mod <- nls(rz~(((1+b*rxn)^(1-c))-1)/(-1+(1+b)^(1-c)), start=list(b=0.001,c=0.1) ,upper = list(b=10,c=10),algorithm = 'port')
    # yhat <- predict(mod,rxn,type="response")
    # plot(rxn,rz)
    # lines(rxn,yhat,col='blue')
    # mod
    # lines(MonoPoly::evalPol(rz,MonoPoly::monpol(rxn~rz,degree = 8,a = 0,b = 1)$beta.raw),rz,col='red')

  }
  return(list('dfnataf'=dfnataf, 'rzEq'=rzEq))
}

TwoParfun<-function(b,c,rx,rmax=1){
  A=(-1+((1+b*rx)^(1-c)))
  B=(-1+(1+b*rmax)^(1-c))
  rzhat=A/B
  return(rzhat)
}

TwoParfunObj<-function(par,rz,rx,rmax=1){
  b=par[1];  c=par[2];
  rzhat=TwoParfun(b = b,c = c,rx = rx,rmax = rmax)

  SSE=sum((rz-rzhat)^2)
  SSE=ifelse(is.nan(SSE),10^6,SSE)
  return(SSE)
}

