#' @title Simulation of correlated Random variables (RVs) with the target marginal distribution and correlation strucutre.
#'
#' @description Simulation of correlated Random variables (RVs) with the target marginal distribution and correlation strucutre.
#'
#' @param n A scalar indicating the number of RVs to generate
#' @param paramsRVs The list that is beeing returned by the function "estCorrRVs"
#'
#' @return A matrix (n x m) with the generated RVs.
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
simCorrRVs=function(n=1000, paramsRVs){
  Rz=paramsRVs$Rz
  b=t(chol(Rz))
  m=ncol(paramsRVs$Rz)
  X = U = matrix(nrow = n, ncol = m)
  V=matrix(data = rnorm(n = n*m), nrow = n, ncol = m)
  Z = t(b %*% t(V))
  
  for (i in 1:m) {
    U[, i] = pnorm(Z[, i])
    fs = paramsRVs$dist[i]
    pfs = paramsRVs$params[[i]]
    X[, i] = eval(as.call(c(as.name(fs), list(U[, i]), pfs)))
  }
  return(X)
}