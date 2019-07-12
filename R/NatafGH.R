#' @title Solve the Nataf integral with the Gauss-Hermite integration method
#'
#' @description Estimation of the resulting (i.e., in the actual domain) correlation coefficients,
#' given the equivalent correlation coefficients (i.e., in the Gaussian domain).
#'
#' @param rho A scalar or vector of correlation coefficients (i.e., in seq(from=0, to=1, by=0.1)).
#' @param fx A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param fy A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param paramlistfx A named list with the parameters of the distribution.
#' @param paramlistfy A named list with parameters of the distribution.
#' @param nodes A scalar indicating the number of nodes for Gauss-Hermite integration (default: nodes=21).
#' @param prune A scalar in (0,1) indicating the percentage of pruning for Gauss-Hermite integration (default: nodes=0).
#'
#' @return A vector of correlation coefficients in the actual domain
#' @note Avoid the use of this function, when the marginal(s) are discrete.
#' @export
#'
#' @examples
#' ## The case of two identrical Gamma distributions, with shape=1 and scale=1.
#'\dontrun{
#' fx=fy='qgamma'
#' pfx=pfy=list(shape=1, scale=1)
#' rhoz=seq(from=0, to=1 , by=0.05)
#' rhox=NatafGH(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx,
#' paramlistfy = pfy, nodes = 10, prune = 0)
#' plot(rhoz,rhox); abline(0,1)
#'
#' ## The case of two identrical zero-inflated (i.e., mixed) distributions,
#' with p0=0.7 a Gamma distribution
#' ## for the continuous part with shape=1 and scale=1.
#'
#' fx=fy='qmixed'
#' pfx=pfy=list(Distr=qgamma, p0=0.7, shape=0.5, scale=1)
#' rhoz=seq(from=0, to=1 , by=0.05)
#' rhox=NatafGH(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx,
#' paramlistfy = pfy, nodes = 21, prune = 0)
#' plot(rhoz,rhox); abline(0,1)
#' }
NatafGH=function(rho, fx, fy, paramlistfx, paramlistfy, nodes=21, prune=0){

  Rx=rho
  Index0=which(rho==0)

  for (t in 1:length(rho)){

    sig <- matrix(data = c(1, rho[t], rho[t], 1), nrow = 2, ncol = 2)
    pts <- mgauss.hermite(nodes, mu=c(0,0), sigma=sig, prune=prune)

    u=pnorm(pts$points)
	  u=ifelse(u==1, 0.999, u)

    funcall1 <- as.call(c(as.name(fx), list(u[, 1]), paramlistfx))
    funcall2 <- as.call(c(as.name(fy), list(u[, 2]), paramlistfy))
    x1=eval(funcall1)
    x2=eval(funcall2)
    Rx[t]=cov.wt(cbind(x1,x2), wt = pts$weights, cor = T, method='ML')$cor[1,2]
    # print(Rx[t])
    Rx[Index0]=0
  }
  return(Rx)
}

## @title Compute the multivariate Gaussian quadrature points
## @description Compute the multivariate Gaussian quadrature points.
##
## @param n Number of points each dimension before pruning.
## @param mu The mean vector.
## @param sigma The covariance matrix.
## @param prune Value specifing the fraction to prune (NULL - no pruning; must be in [0-1])
mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")

  dm  <- length(mu)
  # gh  <- gauss.hermite(n)
  gh=pracma::gaussHermite(n)
  gh$w=gh$w/sqrt(pi) # use this (1)
  gh=cbind(gh$x*sqrt(2),gh$w) # use this (1)
  # gh=cbind(gh$x*sqrt(2),gh$w/sum(gh$w)) # instead of (1) use only this line
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)

  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }

  ## rotate, scale, translate points
  eig <- eigen(sigma)
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}
