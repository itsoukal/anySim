#' @title Simulation of the target cyclostationary process using the SPARTA model of order 1.
#'
#' @description Simulation of the target cyclostationary process using a PAR(p) model (i.e., the SPARTA model of order 1) to simulate the auxiliary cyclostationary Gaussian process to establish the target season-to-season correlation structure.
#'
#' @param SPARTApar A list containing the parameters of the model. The list is constructed by the function "EstSPARTA".
#' @param steps A scalar specifying the length of the time series to be generated.
#' @param stand A boolean (T or F) indicating whether to standardize (or not) the auxiliary Gaussian time series prior to their mapping to the actual domain. The default value is FALSE.
#'
#' @return A list of 3 generated time series (in matrix format - i.e.,  matrix of dimensions k x m, where m denotes the number of sub-seasons and k the number of periods.):
#' X: The final time series at the actual domain with the target marginal distribution and correlation structure;
#' Z: The auxiliary Gaussian time series at the Gaussian domain and,
#' U: The auxiliary uniform time series at the Copula domain (i.e., in [0,1]).
#'
#' @export
#'
#' @examples
#' ## Simulation of univariate cyclostationary process with specific distribution function at each season 
#' ## and specific lag-1 season-to-season correlations.
#' \dontrun{
#' set.seed(21)
#' 
#' # Define the number of seasons.
#' NumOfSeasons=12 # number of months 
#' 
#' # Define the lag-1 season-to-season correlation coefficients (12 values).
#' rtarget<-c(0.05,0.55,0.45,0.4,0.6,0.75,0.7,0.75,0.5,0.3,0.3,0.2)
#' 
#' # Define the target distribution functions for each season.
#' # In this example, a re-parameterized version of Gen. Gamma distribution is used.
#' qgengamma=function(p,scale, shape1, shape2){
#'   require(VGAM)
#'   X=qgengamma.stacy(p=p,scale=scale,k=(shape1/shape2),d=shape2)
#'   return(X)
#' }
#' 
#' # Or, a re-parameterized version of Burr Type XII distribution.
#' qburr=function(p,scale,shape1,shape2) {
#'   require(ExtDist)
#'   x=ExtDist::qBurr(p=p,b=scale,g=shape1,s=shape2)
#'   return(x)
#' }
#' 
#' # Here, we define the target distribution of each season, though being of continuous type,
#' # as a zero-inflated ones to demonstrate the more general case. Alternatively, the definition
#' # could be conducted as in the example of EstARTAp function.
#' 
#' FXs<-rep('qmixed',NumOfSeasons) # Define that distributions are of zero-inflated type.
#' 
#' # Define the parameters of the zero-inflated distribution function for each season.
#' PFXs<-vector("list",NumOfSeasons)
#' PFXs[[1]]=list(p0=0.0,Distr=qgengamma,scale=47.22,shape1=2.7,shape2=0.97)
#' PFXs[[2]]=list(p0=0.0,Distr=qgengamma,scale=199.4,shape1=1.74,shape2=3.45)
#' PFXs[[3]]=list(p0=0.0,Distr=qburr,scale=193.2,shape1=3.07,shape2=2.54)
#' PFXs[[4]]=list(p0=0.0,Distr=qburr,scale=172.16,shape1=4.42,shape2=2.50)
#' PFXs[[5]]=list(p0=0.0,Distr=qgengamma,scale=53.40,shape1=4.11,shape2=1.66)
#' PFXs[[6]]=list(p0=0.0,Distr=qgengamma,scale=0.017,shape1=26.23,shape2=0.51)
#' PFXs[[7]]=list(p0=0.0,Distr=qgengamma,scale=27.70,shape1=5.15,shape2=5.30)
#' PFXs[[8]]=list(p0=0.0,Distr=qgengamma,scale=0.33,shape1=30.97,shape2=0.876)
#' PFXs[[9]]=list(p0=0.0,Distr=qburr,scale=14.46,shape1=7.6,shape2=0.44)
#' PFXs[[10]]=list(p0=0.0,Distr=qburr,scale=29.36,shape1=2.73,shape2=0.87)
#' PFXs[[11]]=list(p0=0.0,Distr=qgengamma,scale=53.15,shape1=3.12,shape2=1.4)
#' PFXs[[12]]=list(p0=0.0,Distr=qgengamma,scale=116.02,shape1=2.21,shape2=1.33)
#'
#' # Estimate the parameters of SPARTA model.
#' SPARTApar<-EstSPARTA(s2srtarget=rtarget,dist=FXs,params=PFXs,
#'             NatafIntMethod='GH',NoEval=9,polydeg=8,nodes=11)
#'             
#' # Generate a synthetic series of 10000 length.
#' simSPARTA<-SimSPARTA(SPARTApar=SPARTApar,steps=10^5)
#' }
SimSPARTA<-function(SPARTApar, steps=1000, stand=0) {
  ## Model Simulation
  NumOfSeasons<-SPARTApar$NumOfSeasons
  FXs=SPARTApar$dist
  PFXs=SPARTApar$params
  r=SPARTApar$s2sEq

  W=matrix(rnorm(steps*NumOfSeasons), steps, NumOfSeasons)
  if (length(r)==1){
    rr=zeros(1,NumOfSeasons)
    rr[1]=r
  }  else {
    rr=r
  }
  r=rr

  m=nrow(W)
  n=ncol(W)

  if (n!=NumOfSeasons) {stop("The number of columns is not equal to NumOfSeasons - The matrix of data is incomplete")}
  Y=matrix(NA,m,n)

  rmod=(1-r^2)^0.5

  for (i in 1:m){
    for (j in 1:n){
      if (j!=1){
        Y[i,j]=Y[i,j-1]*r[j]+rmod[j]*W[i,j]
      } else if (i==1 && j==1) {
        Y[i,j]=W[i,j]
      } else if (j==1) {
        Y[i,j]=Y[i-1,n]*r[j]+rmod[j]*W[i,j]
      }
    }
  }

  X=matrix(NA,m,n)
  U=matrix(NA,m,n)
  for (i in 1:NumOfSeasons) {

    if (stand==0){
      U[,i]=pnorm(Y[,i])
    }else {
      U[,i]=pnorm(Y[,i], mean = mean(Y[,i], sd=sd[Y[,i]]))
    }

    fs=FXs[i]
    pfs=PFXs[[i]]
    X[,i]=eval(as.call(c(as.name(fs), list(U[, i]), pfs)))

  }
  colnames(Y)=paste0('Season_', 1:NumOfSeasons)
  colnames(X)=paste0('Season_', 1:NumOfSeasons)
  colnames(U)=paste0('Season_', 1:NumOfSeasons)
  return(list('Z'=Y, 'X'=X, 'U'=U))
}

