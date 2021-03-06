% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EstSPARTA.R
\name{EstSPARTA}
\alias{EstSPARTA}
\title{Estimation of the auxiliary PAR(1) model parameters}
\usage{
EstSPARTA(
  s2srtarget,
  dist,
  params,
  NatafIntMethod = "GH",
  NoEval = 9,
  polydeg = 6,
  ...
)
}
\arguments{
\item{s2srtarget}{A k-dimensional vector with the lag-1 season-to-season correlation coefficients.}

\item{dist}{A k-dimensional string vector indicating the quantile function of the target marginal distributions (i.e., the ICDFs).}

\item{params}{A k-dimensional named list with the parameters of the target distributions.}

\item{NatafIntMethod}{A string ("GH", "Int", or "MC") indicating the integration method to resolve the Nataf integral.}

\item{NoEval}{A scalar indicating (default: 9) the number of evaluation points for the integration methods.}

\item{polydeg}{A scalar indicating the order of the fitted polynomial. If polydeg=0, then another curve is fitted.}

\item{...}{Additional named arguments for the selected "NatafIntMethod" method.}
}
\value{
A list with the parameters of the cyclostationary auxiliary Gaussian PAR(1) model.
}
\description{
Estimation of parameters of the auxiliary PAR(1) model to simulate the auxiliary cyclostationary Gaussian process to establish the target season-to-season correlation structure.
}
\note{
Avoid the use of the "GH" method (i.e., NatafIntMethod='GH'), when the marginal(s) are discrete.
}
\examples{
## Simulation of univariate cyclostationary process with specific distribution function at each season 
## and specific lag-1 season-to-season correlations.
\dontrun{
set.seed(21)

# Define the number of seasons.
NumOfSeasons=12 # number of months 

# Define the lag-1 season-to-season correlation coefficients (12 values).
rtarget<-c(0.05,0.55,0.45,0.4,0.6,0.75,0.7,0.75,0.5,0.3,0.3,0.2)

# Define the target distribution functions for each season.
# In this example, a re-parameterized version of Gen. Gamma distribution is used.
qgengamma=function(p,scale, shape1, shape2){
  require(VGAM)
  X=qgengamma.stacy(p=p,scale=scale,k=(shape1/shape2),d=shape2)
  return(X)
}

# Or, a re-parameterized version of Burr Type XII distribution.
qburr=function(p,scale,shape1,shape2) {
  require(ExtDist)
  x=ExtDist::qBurr(p=p,b=scale,g=shape1,s=shape2)
  return(x)
}

# Here, we define the target distribution of each season, though being of continuous type,
# as a zero-inflated ones to demonstrate the more general case. Alternatively, the definition
# could be conducted as in the example of EstARTAp function.

FXs<-rep('qmixed',NumOfSeasons) # Define that distributions are of zero-inflated type.

# Define the parameters of the zero-inflated distribution function for each season.
PFXs<-vector("list",NumOfSeasons)
PFXs[[1]]=list(p0=0.0,Distr=qgengamma,scale=47.22,shape1=2.7,shape2=0.97)
PFXs[[2]]=list(p0=0.0,Distr=qgengamma,scale=199.4,shape1=1.74,shape2=3.45)
PFXs[[3]]=list(p0=0.0,Distr=qburr,scale=193.2,shape1=3.07,shape2=2.54)
PFXs[[4]]=list(p0=0.0,Distr=qburr,scale=172.16,shape1=4.42,shape2=2.50)
PFXs[[5]]=list(p0=0.0,Distr=qgengamma,scale=53.40,shape1=4.11,shape2=1.66)
PFXs[[6]]=list(p0=0.0,Distr=qgengamma,scale=0.017,shape1=26.23,shape2=0.51)
PFXs[[7]]=list(p0=0.0,Distr=qgengamma,scale=27.70,shape1=5.15,shape2=5.30)
PFXs[[8]]=list(p0=0.0,Distr=qgengamma,scale=0.33,shape1=30.97,shape2=0.876)
PFXs[[9]]=list(p0=0.0,Distr=qburr,scale=14.46,shape1=7.6,shape2=0.44)
PFXs[[10]]=list(p0=0.0,Distr=qburr,scale=29.36,shape1=2.73,shape2=0.87)
PFXs[[11]]=list(p0=0.0,Distr=qgengamma,scale=53.15,shape1=3.12,shape2=1.4)
PFXs[[12]]=list(p0=0.0,Distr=qgengamma,scale=116.02,shape1=2.21,shape2=1.33)

# Estimate the parameters of SPARTA model.
SPARTApar<-EstSPARTA(s2srtarget=rtarget,dist=FXs,params=PFXs,
            NatafIntMethod='GH',NoEval=9,polydeg=8,nodes=11)
            
# Generate a synthetic series of 10000 length.
simSPARTA<-SimSPARTA(SPARTApar=SPARTApar,steps=10^5)
}
}
