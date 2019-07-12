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
#' ## Simulation example of a bivariate process with zero-inflated marginal distributions.
#' # Define the simulation parameters ----------------------------------------
#' \dontrun{
#' LAG=2^6
#' FFTLag=2^7
#' SMALAG=2^6
#' steps=2^14
#'
#' PFXs=list()
#' FXs=c('qmixed','qmixed')
#' # Gamma distribution: Gamma(shape=2, rate=1)
#' PFXs[[1]]=list(Distr=qgamma, p0=0.9, shape=1, scale=1)
#' # Weibull distribution: Weibull(shape=1, scale=2)
#' PFXs[[2]]=list(Distr=qweibull, p0=0.85, shape=1, scale=2)
#'
#' ACFs=list()
#' ACFs[[1]]=acsCAS(param = c(0.1, 0.6), lag = LAG)
#' ACFs[[2]]=acsCAS(param = c(0.2, 0.3), lag = LAG)
#'
#' Cmat=matrix(c(1,0.6,0.6,1), ncol=2, nrow=2)
#'
#' # Calculate SMARTA's parameters -------------------------------------------
#' SMAparam=EstSMARTA(dist = FXs, params = PFXs, ACFs = ACFs, Cmat = Cmat,
#' DecoMethod = 'cor.smooth',FFTLag = FFTLag,
#' NatafIntMethod = 'GH', NoEval = 9, polydeg = 8)
#'
#' # Simulate a SMARTA process -----------------------------------------------
#' Sim=SimSMARTA(SMARTApar = SMAparam, steps = steps, SMALAG = SMALAG)
#'
#' # Draw some basic plots ---------------------------------------------------
#' for (i in 1:2) {
#'  par(mfrow=c(2,2))
#'   plot(Sim$X[1:1001,i], type='l')
#'  acf(Sim$X[,i], lag.max = 20); lines(0:20,ACFs[[i]][1:21], col='red', type='o')
#'  plot(ecdf(Sim$X[,i]))
#'   hist(Sim$X[,i],probability = TRUE)
#' }
#'}
#'
#'\dontrun{
#' ## Simulation example of a bivariate process with Bernoulli marginal distributions.
#' # Define the simulation parameters ----------------------------------------
#' PFXs=list()
#' FXs=c('qbinom','qbinom')
#' PFXs[[1]]=list(size=1, prob=0.2)# Gamma distribution: Gamma(shape=2, rate=1)
#' PFXs[[2]]=list(size=1, prob=0.25) # Weibull distribution: Weibull(shape=1, scale=2)
#'
#' ACFs=list()
#' ACFs[[1]]=acsCAS(param = c(0.1, 0.6), lag = LAG)
#' ACFs[[2]]=acsCAS(param = c(0.2, 0.3), lag = LAG)
#'
#' Cmat=matrix(c(1,0.6,0.6,1), ncol=2, nrow=2)
#'
#' # Calculate SMARTA’s parameters -------------------------------------------
#' SMAparam=EstSMARTA(dist = FXs, params = PFXs, ACFs = ACFs, Cmat = Cmat, DecoMethod = 'cor.smooth',
#'                    FFTLag = FFTLag, NatafIntMethod = 'Int', NoEval = 9, polydeg = 8)
#' # Simulate a SMARTA process -----------------------------------------------
#' Sim=SimSMARTA(SMARTApar = SMAparam, steps = steps, SMALAG = SMALAG)
#'
#' # Draw some basic plots ---------------------------------------------------
#' for (i in 1:2) {
#'   par(mfrow=c(2,2))
#'   plot(Sim$X[1:1001,i], type='l')
#'   acf(Sim$X[,i], lag.max = 20); lines(0:20,ACFs[[i]][1:21], col='red', type='o')
#'   plot(ecdf(Sim$X[,i]))
#'   hist(Sim$X[,i],probability = TRUE)
#' }
#'}
#'
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
