#' @title Simulation of the target stationary processes using the SMARTA model.
#'
#' @description Simulation of the target stationary processes using an SMARTA model to simulate the auxiliary Gaussian process and establish the target correlation structure.
#' This model allows also the simulation of a multivariate process.
#'
#' @param SMARTApar A list containing the parameters of the model. The list is constructed by the function "EstSMARTA".
#' @param steps A scalar specifying the length of the time series to be generated.
#' @param SMALAG A scalar specifying the order of the SMARTA model (must be equal to FFTAG and the lengths of ACFs).
#'
#' @return A list of the 3 simulated time series (in matrix format - i.e.,  matrix of dimensions steps x k, where steps denotes the length of generated time series and k the number of processes).
#' X: The final time series at the actual domain with the target marginal distribution and correlation structure;
#' Z: The auxiliary Gaussian time series at the Gaussian domain and;
#' U: The auxiliary uniform time series at the Copula domain (i.e., in [0,1]).
#'
#' @export
#'
#' @examples
#' ## Multivariate simulation of 3 stationary processes with specific distribution functions 
#' ## and autocorrelation structures, as well as specific lag-0 cross-correlation matrix.
#' \dontrun{
#' set.seed(9)
#' 
#' # Define the target autocorrelation structure of the 3 processes.
#' ACSs=list()
#' ACSs[[1]]=csCAS(param=c(0.1,0.7),lag=2^6)
#' ACSs[[2]]=csCAS(param=c(0.2,1),lag=2^6)
#' ACSs[[3]]=csCAS(param=c(0.1,0.5),lag=2^6)
#' 
#' # Define the matrix of lag-0 cross-correlation coefficients between the 3 processes.
#' Cmat=matrix(c(1,0.4,-0.5,
#'               0.4,1,-0.3,
#'              -0.5,-0.3,1),ncol=3,nrow=3)
#' 
#' # Define the target distribution functions (ICDF) of the 3 processes
#' FXs=rep('qmixed',3) # Define that distributions are of zero-inflated type.
#' 
#' # Define the distributions for the continuous part of the processes. 
#' # In this example, a re-parameterized version of Gen. Gamma distribution is used for the second process.
#' qgengamma=function(p,scale, shape1, shape2){
#'   require(VGAM)
#'   X=qgengamma.stacy(p=p,scale=scale,k=(shape1/shape2),d=shape2)
#'   return(X)
#' }
#' 
#' # Define the parameters of the target distributions.
#' pFXs[[1]]=list(Distr=qbeta,p0=0,shape1=15,shape2=5) # Beta distribution
#' pFXs[[2]]=list(Distr=qgengamma,p0=0.7,scale=0.12, shape1=1.35, shape2=0.4) # Gen. Gamma
#' pFXs[[3]]=list(Distr=qnorm,p0=0,mean=15,sd=3) # Normal distribution
#' 
#' # Estimate the parameters of SMARTA model
#' SMAparam=EstSMARTA(dist=FXs,params=pFXs,ACFs=ACSs,Cmat=Cmat,
#'                    DecoMethod='cor.smooth',FFTLag = 2^7,
#'                    NatafIntMethod='GH',NoEval=9,polydeg=8)
#'                    
#' # Generate the synthetic series     
#' simSMARTA=SimSMARTA(SMARTApar=SMAparam,steps=2^14,SMALAG=2^6)
#'}
SimSMARTA <- function(SMARTApar, steps, SMALAG=512) {


  YEARSTOT=2*SMALAG+1+steps
  Sites=SMARTApar$Sites

  # define list of synthetic data for each site
  SMARTApar.names=paste("Site", 1:Sites,sep="")
  Xlist=vector("list", length(SMARTApar.names))

  X=U=Z=matrix(nrow=steps, ncol=Sites) # matrix of synthetic data
  V=matrix(nrow=YEARSTOT, ncol=Sites)

  # Generate uncorrelated Innovations V
  for (i in 1:Sites)  {
    V[,i]=rnorm(n = YEARSTOT)
  }

  # Produce cross-correlated Innovations W
  W=t(SMARTApar$b%*%t(V))

  for (i in 1:Sites) {
    aj=t(SMARTApar$FFTa[[i]])
    X[,i]=sumprod(aj=aj, W=as.matrix(W[,i], ncol=1), Years=steps, SMALAG=SMALAG)
  }
  Z=X
  # Map the auxiliary data to the actual domain
  for (i in 1:Sites) {
    U[,i]=pnorm(X[,i])
    fs=SMARTApar$icdfname[i]
    pfs=SMARTApar$params[[i]]
    # X[,i]=do.call(fs, args = c(U,pfs))
    X[,i]=eval(as.call(c(as.name(fs), list(U[,i]), pfs)))
  }

  # assign synthetic data for each site to the list
  for (i in 1:Sites) {
    Xlist[[i]]=matrix(X[,i], ncol=1)
  }

  return(list('Z'=Z, 'X'=X, 'U'=U))
}

sumprod <- function(aj, W, Years, SMALAG) {

  # define matrix of synthetic data
  X=matrix(nrow=Years,ncol=1)

  # the sumproducts of SMA scheme
  for (i in 1:Years) {
    endindex=2*SMALAG+i
    MAT=W[i:endindex,1]
    X[i,1]=sum(aj*MAT)
  }

  X=as.vector(X)
  return(X)

}
