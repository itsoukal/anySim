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
#' ## Simulation of cyclostationary process with 12 seasons and zero-inflated marginal distributions.
#' \dontrun{
#' rtarget<-c(0.5, 0.7, 0.6, 0.4, 0.5, 0.7, 0.8, 0.7, 0.6, 0.4, 0.5, 0.7)
#' NumOfSeasons=length(rtarget)
#' FXs<-rep('qmixed',NumOfSeasons)
#' PFXs<-vector("list",NumOfSeasons)
#' PFXs<-lapply(PFXs,function(x) x<-list(p0=0.4, Distr=qexp, rate=0.5))
#'
#' SPARTApar<-EstSPARTA(s2srtarget = rtarget, dist = FXs, params = PFXs,
#' NatafIntMethod = 'GH', NoEval = 9, polydeg = 8)
#'
#' sim<-SimSPARTA(SPARTApar = SPARTApar, steps=100000, stand=0)
#' s2scor(sim$X)
#' plot(s2scor(sim$X)); lines(rtarget, col='red')
#'}
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

