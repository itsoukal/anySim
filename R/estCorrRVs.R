#' @title Estimation of the auxiliary Gaussian model parameters for the generation of correlated random variables.
#'
#' @description Estimation of the auxiliary Gaussian model parameters for the generation of correlated Random variables (RVs).
#'
#' @param R A k x k matrix with the target Pearson correlation coefficients among the k RVs.
#' @param dist A k-dimensional string vector indicating the quantile function of the target marginal distributions (i.e., the ICDF) of k RVs.
#' @param params A k-dimensional named list with the parameters of the k target distributions.
#' @param NatafIntMethod A string ("GH", "Int", or "MC") indicating the intergation method to resolve the Nataf integral.
#' @param NoEval A scalar indicating (default: 9) the number of evaluation points for the integration methods.
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
#' ## Simulation of 3 correlated RVs with Gamma, Beta and Log-Normal distribution, respectively. 
#' ## We assume also a target correlation matrix.
#' \dontrun{
#' set.seed(13)
#' 
#' # Define the target distribution functions (ICDFs) of X1, X2 and X3 RV.
#' FX1='qgamma'; FX2='qbeta'; FX3='qlnorm'
#' Distr=c(FX1,FX2,FX3) # store in a vector.
#'
#' # Define the parameters of the target distribution functions 
#' # and store them in a list.  
#' pFX1=list(shape=1.5,scale=2); pFX2=list(shape1=1.5,shape2=3)
#' pFX3=list(meanlog=1,sdlog=0.5)
#' 
#' DistrParams=list()
#' DistrParams[[1]]=pFX1;DistrParams[[2]]=pFX2;DistrParams[[3]]=pFX3
#' 
#' # Define the target correlation matrix.
#' CorrelMat=matrix(c(1,0.7,0.5,
#'                    0.7,1,0.8,
#'                    0.5,0.8,1),ncol=3,nrow=3,byrow=T);
#'             
#' # Estimate the parameters of the auxiliary Gaussian model.
#' paramsRVs=EstCorrRVs(R=CorrelMat,dist=Distr,params=DistrParams, 
#'                         NatafIntMethod='GH',NoEval=9,polydeg=8)
#'
#' # Generate 10000 synthetic realisations of the 3 correlated RVs.
#' SynthRVs=SimCorrRVs(n=10000,paramsRVs=paramsRVs)
#' 
#' # Comparison of target and simulation correlation matrix.
#' pairs(SynthRVs)
#' CorrelMatSim=cor(SynthRVs)
#' difference=(CorrelMat-CorrelMatSim)^2; print(difference)
#'}
EstCorrRVs=function(R, dist, params, NatafIntMethod = 'GH', NoEval = 9, polydeg = 8) {
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
