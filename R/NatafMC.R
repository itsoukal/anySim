#' @title Solve the Nataf integral with the Monte-Carlo method
#'
#' @description Estimation of the resulting (i.e., in the actual domain) correlation coefficients,
#' given the equivalent correlation coefficients (i.e., in the Gaussian domain).
#'
#' @param rho A scalar or vector of correlation coefficients (i.e., in seq(from=0, to=1, by=0.1)).
#' @param fx A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param fy A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param paramlistfx A named list with the parameters of the distribution.
#' @param paramlistfy A named list with parameters of the distribution.
#' @param MCsize A scalar determining the number of Monte-Carlo trials (default: 10^-5)
#' @return A vector of correlation coefficients in the actual domain.
#' @export
#'
#' @examples
#' ## The case of two identical Gamma distributions, with shape=1 and scale=1.
#'\dontrun{
#' fx=fy='qgamma'
#' pfx=pfy=list(shape=1, scale=1)
#' rhoz=seq(from=0, to=1 , by=0.2)
#' rhoxMC=NatafMC(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' plot(rhoz,rhoxMC, col='blue', pch=19); abline(0,1)
#' rhoxGH=NatafGH(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' points(rhoz, rhoxGH, col='red', pch=5)
#'
#' ## The case identrical Bernoulli distributions, with size=1 and prob=0.3.
#'
#' fx=fy='qbinom'
#' pfx=pfy=list(size=1, prob=0.3)
#' rhoz=seq(from=0, to=1 , by=0.2)
#' rhoxMC=NatafMC(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' plot(rhoz,rhoxMC, col='blue', pch=19); abline(0,1)
#'}
#'
#' \dontrun{
#' rhoxInt=NatafInt(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' points(rhoz, rhoxInt, col='red', pch=5)
#'}
NatafMC=function(rho, fx, fy, paramlistfx, paramlistfy, MCsize=10^5){

  Rx=rho
  Index0=which(rho==0)

  # Z=MASS::mvrnorm(n = MCsize, mu = c(0,0), Sigma = matrix(c(1,0,0,1),2,2))
  Z=cbind(rnorm(MCsize),rnorm(MCsize))
  ZZ=matrix(NA, nrow=MCsize, ncol=2)
  for (t in 1:length(rho)){
    ZZ[,1]=Z[,1]
    ZZ[,2]=Z[,1]*rho[t]+Z[,2]*sqrt(1-rho[t]^2)
    U=pnorm(ZZ)
    funcall1 <- as.call(c(as.name(fx), list(U[, 1]), paramlistfx))
    funcall2 <- as.call(c(as.name(fy), list(U[, 2]), paramlistfy))
    x1=eval(funcall1)
    x2=eval(funcall2)
    Rx[t]=cor(x1,x2)
  }
  Rx[Index0]=0
  return(Rx)
}

