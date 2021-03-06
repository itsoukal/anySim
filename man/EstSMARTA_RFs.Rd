% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EstSMARTA_RFs.R
\name{EstSMARTA_RFs}
\alias{EstSMARTA_RFs}
\title{Estimation of the auxiliary SMA model parameters}
\usage{
EstSMARTA_RFs(
  dist,
  params,
  ACFs,
  Cmat,
  DecoMethod = "cor.smooth",
  FFTLag = 512,
  NatafIntMethod = "GH",
  NoEval = 9,
  polydeg = 8,
  ...
)
}
\arguments{
\item{dist}{A k-dimensional string vector indicating the quantile function of the target marginal distribution (i.e., the ICDF).}

\item{params}{A k-dimensional named list with the parameters of the target distributions.}

\item{ACFs}{A k-dimensional list with the target autocorrelation structure (including lag-0, i.e., 1).}

\item{Cmat}{A matrix (k x k) containing the lag-0 cross-correlation coefficients of the processes.}

\item{DecoMethod}{A string indicating the decomposition method, in case of a non-positive definite matrix (options: 'cor.smooth' and 'nearPD')}

\item{FFTLag}{A scalar indicating the length of the Fast Fourrier Transform (required to estimate the internal parameters of SMA model). Default value=512.}

\item{NatafIntMethod}{A string ("GH", "Int", or "MC"), indicating the intergation method, to resolve the Nataf integral.}

\item{NoEval}{A scalar (power of 2) indicating (default: 9) the number of evaluation points for the integration methods.}

\item{polydeg}{A scalar indicating the order of the fitted polynomial. If polydeg=0, then another curve is fitted.}

\item{...}{Additional named arguments for the selected "NatafIntMethod" method.}
}
\value{
A list with the parameters of the auxiliary Gaussian SMA model.
}
\description{
Estimation of parameters of the auxiliary SMA model to simulate the auxiliary Gaussian process.
}
\note{
Avoid the use of the "GH" method (i.e., NatafIntMethod='GH'), when the marginal(s) are discrete.
}
\examples{
## Simulation of homogenous and stationary non-Gaussian random fields (RFs)
\dontrun{

# Define a 30x30 grid to be simulated
nx=30 # number of cells in the horizontal direction
ny=30 # number of cells in the vertical direction
Sites=nx*ny # number of grid points

Xp=seq(from=(0.5),to=nx,by=1) # points' coordinates in horizontal axis
Yp=seq(from=(0.5),to=ny,by=1) # points' coordinates in vertical axis

grid=expand.grid(X=Xp,Y=Yp)

# plot(grid,cex=1,pch=19,col='lightblue',xlim=c(0, nx),ylim=c(0, ny))
# text(grid,labels=rownames(grid),cex=0.5,font=2)
# abline(h=0:nx,v=0:ny)
# plot(grid,pch=19)

# Estimate the Euclidean distances between grid points
DZ=dist(x=grid,method='euclidean',upper=T,diag=T)
DZmat=as.matrix(DZ)
EuclDist=DZmat[upper.tri(DZmat, diag = T)]

# Define the matrix of lag-0 cross-correlation coefficients among grid points.
CCF=(1+0.2*2*EuclDist)^(-1/b) # CAS with b=0.2 and k=2.
Cmat=matrix(NA,nrow=nx*ny,ncol=nx*ny)
Cmat[upper.tri(Cmat,diag=T)]=CCF
Cmat[lower.tri(Cmat,diag=T)]=rev(CCF)

# Define the target autocorrelation structure and distribution function (ICDF) at each point.
# The distribution is of zero-inflated type with Burr Type-XII distribution for the continuous part.
# Here, the following re-parameterized version of Burr Type-XII distribution is used.
qburr=function(p,scale,shape1,shape2) {
      require(ExtDist)
      x=ExtDist::qBurr(p=p,b=scale,g=shape1,s=shape2)
      return(x)
}

FXs=rep('qzi',Sites) # Define that distributions are of zero-inflated type.
PFXs=vector("list",length=Sites) # List with ICDF of each point
ACFs=vector("list",length=Sites) # List with ACF of each point
for (i in 1:Sites) {
 PFXs[[i]]=list(Distr=qburr,p0=0.75,scale=71.62,shape1=0.88,shape2=11.79)
 ACFs[[i]]=csCAS(param=c(0.1,0.6),lag=2^6) # CAS with b=0.1 and k=0.6.
}

# Estimate the parameters of SMARTA model
SMAparam=EstSMARTA_RFs(dist=FXs,params=PFXs,ACFs=ACFs,
                       Cmat=Cmat,DecoMethod='cor.smooth',
                       FFTLag=2^7,NatafIntMethod='GH',NoEval=9,polydeg=8)

# Generate a synthetic realisation of random fields with 2^15 length     
SimField=SimSMARTA(SMARTApar=SMAparam,steps=2^15,SMALAG=2^6)
}
}
