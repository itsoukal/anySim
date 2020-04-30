#' @title A classical periodic autocorrelation structure
#' @description A classical periodid autocorrelation structure derived by David Mackay.
#' @param param A two-dimenstional vector, containing the parameters of the periodic autocorrelation structure. First position is period (p), and the second the lengthscale (l).
#' @param lag A scalar indicating the maximum lag, up to which the fGn is estimated.
#' @param var A scalar indicating the variance of the process. If var!=1 then, the autocovariance structure is returned, instead of the autocorrelation structure.
#'
#' @return A vector of length (lag+1) with the values of the periodic autocorrelation structure.
#' @export
#'
#' @examples
#' ## Estimate and plot a periodic autocorrelation structure with parameters, period=12 and
#' lenghtscale=1, up to lag 50.
#'
#' ACF=acsPeriodic(param=c(12, 1), lag=50, var=1)
#' plot(0:(length(ACF)-1), ACF, t='l')
csPeriodic=function(param, lag, var=1) {
  d=0:lag
  A=pi*d/param[1]
  rho=exp(-2*sin(A)^2/param[2]^2)*var
  return(rho)
}
