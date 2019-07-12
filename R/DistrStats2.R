#' @title Estimation of the theoretical Mean, Variance, Skewness and Kurtosis coefficients of a distribution function (from the quantile function).
#'
#' @description Estimate the theoretical Mean, Variance, Skewness and Kurtosis coefficients of a distribution function, based on the quantile function.
#'
#' @param qDistr The quantile function (as a function) of the distribution.
#' @param lb A scalar indicating the lower bound of distribution.
#' @param ub A scalar indicating the upper bound of distribution.
#' @param subdiv A scalar indicating the maximum number of subintervals.
#' @param rel.tol A scalar indicating the relative accuracy requested.
#' @param abs.tol A scalar indicating the absolute accuracy requested.
#' @param ... Additional named arguments containing the distribution parameters.
#'
#'
#' @return A 4-dimensional vector containing the theoretical Mean, Variance, Skewness and Kurtosis coefficients of the distribution.
#' @export
#'
#' @examples
#' ## Gamma distribution with shape=0.5 and scale=2
#' DistrStats2(qDistr = qgamma, shape=1,scale=1, lb = 0,ub = 1)
#'
#' ## Zero-inflated (i.e., mixed) distribution with p0=0.7 and
#' ## continuous part given by Gamma distribution with shape=0.5 and scale=2
#'
#' DistrStats2(qDistr = function(x) qmixed(x,Distr = qgamma, p0=0.7, shape=1, scale=1),lb = 0,ub=1)
#'
DistrStats2<-function(qDistr, lb=0, ub=1, subdiv=90000, rel.tol = 10^-8, abs.tol=10^-8,  ...){

  mu1=integrate(f=(function(x) qDistr(x,...)),lower = lb,upper = ub, subdivisions = subdiv, rel.tol = rel.tol, abs.tol = abs.tol)$value # mean
  var=(integrate(f=(function(x) qDistr(x,...)^2),lower = lb,upper = ub, subdivisions = subdiv, rel.tol = rel.tol, abs.tol = abs.tol)$value-mu1^2) # var
  sigma1=sqrt(integrate(f=(function(x) qDistr(x,...)^2),lower = lb,upper = ub,subdivisions = subdiv, rel.tol = rel.tol, abs.tol = abs.tol)$value-mu1^2) # sd
  sk1=(integrate(f=(function(x) (qDistr(x,...)-mu1 )^3 ),lower = lb,upper = ub, subdivisions = subdiv, rel.tol = rel.tol, abs.tol = abs.tol)$value)/sigma1^3 # skewness
  kurt1=(integrate(f=(function(x) (qDistr(x,...)-mu1)^4),lower = lb,upper = ub,subdivisions = subdiv, rel.tol = rel.tol, abs.tol = abs.tol)$value)/sigma1^4 # not excess kurtosis

  stats<-c(mu1,var,sk1,kurt1)
  names(stats)<-c("Mean","Var","Skewness","Kurtosis")

  return(stats)

}
