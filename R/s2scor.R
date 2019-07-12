#' @title Estimation of lag-1 season-to-season correlation coefficients
#'
#' @description Estimation of lag-1 seasson-to-season correlation coefficients.
#'
#' @param data A matrix of dimensions k x m, where m denotes the number of sub-seasons and k the number of periods.
#' For instance k may refer to years and m to months (i.e., 12). Another example could regard k as days (e.g., 31 x years) and m to hours (i.e., 24).
#'
#' @return A k-dimensional vector with the lag-1 season-to-season correlations coefficients.
#' @export
#'
#' @examples
#' ## Simulation of cyclostationary process with 12 seasons and zero-inflated marginal distributions.
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
#'
s2scor<-function(data) {
  TempCor=cor(data,use='pairwise.complete.obs')
  v1=cor(data[-1,1],data[-nrow(data),ncol(data)],use='pairwise.complete.obs')
  v2=diag(TempCor[-1,])

  corLag1=round(c(v1,v2),6)
  return(corLag1)
}
