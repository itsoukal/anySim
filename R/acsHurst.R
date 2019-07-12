#' @title The Hurst or fGn autocorrelation structure
#' @description The Hurst or fGn autocorrelation structure
#' @param H A scalar indicating the Hurst coefficient.
#' @param lag A scalar indicating the maximum lag, up to which the fGn is estimated.
#' @param var A scalar indicating the variance of the process. If var!=1 then, the autocovariance structure is returned, instead of the autocorrelation structure.
#'
#' @return A vector of length (lag+1) with the values of fGn autocorrelation structure.
#' @export
#'
#' @examples
#' ## Estimate and plot an fGn (i.e., Hurst) structure with H=0.7, up to lag 500.
#' ACF=acsHurst(H=0.7, lag=500, var=1)
#' plot(0:(length(ACF)-1), ACF)
acsHurst<- function(H,lag,var=1){
  g=rep(NA,lag)
  g1=1

  for (i in 1:lag){
    g[i]=(1/2)*((i+1)^(2*H)+(i-1)^(2*H))-i^(2*H)

  }
  g=(c(g1,g))*var
  return(g)

}
