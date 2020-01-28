#' @title Estimation of the auxiliary Gaussian model parameters for the generation of correlated Random variables (RVs).
#'
#' @description Estimation of the auxiliary Gaussian model parameters for the generation of correlated Random variables (RVs).
#'
#' @param R A k x k matrix with the target pearson correlation coefficients.
#' @param dist A k-dimensional string vector indicating the quantile function of the target marginal distribution (i.e., the ICDF).
#' @param params A k-dimensional named list with the parameters of the target distributions.
#' @param NoEval A scalar (power of 2) indicating (default: 9) the number of evaluation points for the integration methods.
#' @param polydeg A scalar indicating the order of the fitted polynomial. If polydeg=0, then another curve is fitted.
#' @param ... Additional named arguments for the selected "NatafIntMethod" method.
#'
#' @note Avoid the use of the "GH" method (i.e., NatafIntMethod='GH'), when the marginal(s) are discrete.
#'
#' @return A list with the parameters of the auxiliary Gaussian  model.
#'
#' @export
#'
#' @examples
#' ## Simulation example of a bivariate process with
#' ## zero-inflated marginal distributions.
#' # Define the simulation parameters ----------------------------------------
#' \dontrun{
#' Define the target distribution function for X1, X2, and X3.
#' FX1='qgamma'
#' FX2='qbeta'
#' FX3='qlnorm'
#' dist=c(FX1, FX2, FX3) # Store them in a vector.
#' 
#' # Define the parameters of the target distribution functions.
#' pFX1=list(shape=1.5, scale=2)
#' pFX2=list(shape1=1.5, shape2=3)
#' pFX3=list(meanlog=1, sdlog=0.5)
#' params=list() # Store them in a vector.
#' params[[1]]=pFX1
#' params[[2]]=pFX2
#' params[[3]]=pFX3
#' 
#' # Define the target correlation matrix.
#' R=matrix(data = c(1, 0.7, 0.5, 
#'                   0.7, 1, 0.8, 
#'                   0.5, 0.8, 1), ncol=3, nrow = 3, byrow = T);R
#' paramsRVs=estCorrRVs(R = R, dist = dist, params = params, NatafIntMethod = 'GH', NoEval = 9, polydeg = 8)
#' 
#' X=simCorrRVs(n = 10000, paramsRVs = paramsRVs)
#' pairs(X)
#' Rsim=cor(X)
#' difference=(R-Rsim)^2; print(difference)
#'}
estCorrRVs=function(R, dist, params, NatafIntMethod = 'GH', NoEval = 9, polydeg = 8) {
  m=ncol(R)
  Rz=matrix(data = NA, nrow = m, ncol = m); diag(Rz)=1
  for (i in 1:(m-1)) {
    for (j in (i+1):m){
      Rz[i,j]=Rz[j,i]=NatafInvD(targetrho = R[i, j], fx = dist[i], fy = dist[j], 
                                paramlistfx = params[[i]], paramlistfy = params[[j]], 
                                NatafIntMethod = NatafIntMethod, NoEval = NoEval, polydeg = polydeg)$rzEq
    }
  }
  return(list('Rz'=Rz, 'dist'=dist, 'params'=params))
}