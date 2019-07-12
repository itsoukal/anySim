#' @title The Cauchy-type autocorrelation structure (CAS)
#' @description The Cauchy-type autocorrelation structure (CAS)
#' @param param A two-dimenstional vector, containing the parameters of CAS. First position is for b, the second for k.
#' @param lag A scalar indicating the maximum lag, up to which CAS is estimated.
#' @param var A scalar indicating the variance of the process. If var!=1 then,the autocovariance structure is returned, instead of the autocorrelation structure.
#'
#' @return A vector of length (lag+1) with the values of CAS autocorrelation structure.
#' @export
#'
#' @examples
#' ##  Estimate and plot a CAS structure, with b=1 and k=0.6, up to lag 500.
#' ACF=acsCAS(param=c(1, 0.6), lag=500, var=1)
#' plot(0:(length(ACF)-1), ACF)
acsCAS <- function(param,lag,var=1) {

  b=param[1];k=param[2]

  i=1:lag
  if(b==0){
    correl=exp(-k*i)
  } else {
    correl=(1+k*b*i)^(-1/b)
  }

  correl=c(1,correl)
  correl=as.vector(var*t(correl))

  return(correl)
}
