#' @title Estimation of the auxiliary SMA model parameters
#'
#' @description Estimation of parameters of the auxiliary SMA model to simulate the auxiliary Gaussian process.
#'
#' @param dist A k-dimensional string vector indicating the quantile function of the target marginal distribution (i.e., the ICDF).
#' @param params A k-dimensional named list with the parameters of the target distributions.
#' @param ACFs A k-dimensional list with the target autocorrelation structure (including lag-0, i.e., 1).
#' @param Cmat A matrix (k x k) containing the lag-0 cross-correlation coefficients of the processes.
#' @param DecoMethod A string indicating the decomposition method, in case of a non-positive definite matrix (options: 'cor.smooth' and 'nearPD')
#' @param FFTLag A scalar indicating the length of the Fast Fourrier Transform (required to estimate the internal parameters of SMA model). Default value=512.
#' @param NatafIntMethod A string ("GH", "Int", or "MC"), indicating the intergation method, to resolve the Nataf integral.
#' @param NoEval A scalar (power of 2) indicating (default: 9) the number of evaluation points for the integration methods.
#' @param polydeg A scalar indicating the order of the fitted polynomial. If polydeg=0, then another curve is fitted.
#' @param ... Additional named arguments for the selected "NatafIntMethod" method.
#'
#' @note Avoid the use of the "GH" method (i.e., NatafIntMethod='GH'), when the marginal(s) are discrete.
#'
#' @return A list with the parameters of the auxiliary Gaussian SMA model.
#'
#' @export
#'
#' @examples
#' ## Simulation of homogenous and stationary non-Gaussian random fields (RFs)
#' \dontrun{
#' 
#' # Define a 30x30 grid to be simulated
#' nx=30 # number of cells in the horizontal direction
#' ny=30 # number of cells in the vertical direction
#' Sites=nx*ny # number of grid points
#' 
#' Xp=seq(from=(0.5),to=nx,by=1) # points' coordinates in horizontal axis
#' Yp=seq(from=(0.5),to=ny,by=1) # points' coordinates in vertical axis
#' 
#' grid=expand.grid(X=Xp,Y=Yp)
#' 
#' # plot(grid,cex=1,pch=19,col='lightblue',xlim=c(0, nx),ylim=c(0, ny))
#' # text(grid,labels=rownames(grid),cex=0.5,font=2)
#' # abline(h=0:nx,v=0:ny)
#' # plot(grid,pch=19)
#' 
#' # Estimate the Euclidean distances between grid points
#' DZ=dist(x=grid,method='euclidean',upper=T,diag=T)
#' DZmat=as.matrix(DZ)
#' EuclDist=DZmat[upper.tri(DZmat, diag = T)]
#' 
#' # Define the matrix of lag-0 cross-correlation coefficients among grid points.
#' CCF=(1+0.2*2*EuclDist)^(-1/b) # CAS with b=0.2 and k=2.
#' Cmat=matrix(NA,nrow=nx*ny,ncol=nx*ny)
#' Cmat[upper.tri(Cmat,diag=T)]=CCF
#' Cmat[lower.tri(Cmat,diag=T)]=rev(CCF)
#' 
#' # Define the target autocorrelation structure and distribution function (ICDF) at each point.
#' # The distribution is of zero-inflated type with Burr Type-XII distribution for the continuous part.
#'# Here, the following re-parameterized version of Burr Type-XII distribution is used.
#' qburr=function(p,scale,shape1,shape2) {
#'       require(ExtDist)
#'       x=ExtDist::qBurr(p=p,b=scale,g=shape1,s=shape2)
#'       return(x)
#' }
#' 
#' FXs=rep('qzi',Sites) # Define that distributions are of zero-inflated type.
#' PFXs=vector("list",length=Sites) # List with ICDF of each point
#' ACFs=vector("list",length=Sites) # List with ACF of each point
#' for (i in 1:Sites) {
#'  PFXs[[i]]=list(Distr=qburr,p0=0.75,scale=71.62,shape1=0.88,shape2=11.79)
#'  ACFs[[i]]=csCAS(param=c(0.1,0.6),lag=2^6) # CAS with b=0.1 and k=0.6.
#' }
#' 
#' # Estimate the parameters of SMARTA model
#' SMAparam=EstSMARTA_RFs(dist=FXs,params=PFXs,ACFs=ACFs,
#'                        Cmat=Cmat,DecoMethod='cor.smooth',
#'                        FFTLag=2^7,NatafIntMethod='GH',NoEval=9,polydeg=8)
#' 
#' # Generate a synthetic realisation of random fields with 2^15 length     
#' SimField=SimSMARTA(SMARTApar=SMAparam,steps=2^15,SMALAG=2^6)
#'}
EstSMARTA_RFs<-function(dist, params, ACFs, Cmat, DecoMethod="cor.smooth", FFTLag=512, NatafIntMethod='GH', NoEval=9, polydeg=8, ...) {

  Sites=length(dist)
  # GAFLag=512     # Number of lags used to calc theoretical GAF
  # FFTLag=512      # Number of lags used in FFT procedure

  SMAparams=list()
  SMAparams$Sites=Sites
  # Define arguments of "param" list (sub-lists)
  namesSites=paste("Site",1:Sites,sep="_")
  NatafdflistAutocorr=list()

  for (i in 1:1) {
    # store distribution parameters
    SMAparams$icdfname=dist
    SMAparams$params[[i]]=params[[i]]

    paramslist=list()
    paramslist$m1=params[[i]]
    paramslist$m2=params[[i]]

    if (NatafIntMethod=='GH') {
      Nataf=NatafInvD(targetrho = ACFs[[i]][-1], fx = dist[i], fy = dist[i], paramlistfx = params[[i]], paramlistfy = params[[i]],
                      NatafIntMethod = 'GH', polydeg=polydeg, NoEval = NoEval, ...)
    } else if (NatafIntMethod=='Int') {
      Nataf=NatafInvD(targetrho = ACFs[[i]][-1], fx = dist[i], fy = dist[i], paramlistfx = params[[i]], paramlistfy = params[[i]],
                      NatafIntMethod = 'Int', polydeg=polydeg, NoEval = NoEval)
    } else if (NatafIntMethod=='MC') {
      Nataf=NatafInvD(targetrho = ACFs[[i]][-1], fx = dist[i], fy = dist[i], paramlistfx = params[[i]], paramlistfy = params[[i]],
                      NatafIntMethod = 'MC', polydeg=polydeg, NoEval = NoEval, ...)
    } else {
      print('Select an appropriate method to resolve Nataf integral. That is: "GH, "MC" or "Int"')
    }

    NatafdflistAutocorr[[i]]=Nataf$dfnataf

    SMAparams$GAFv[[i]]=(c(1,Nataf$rzEq))
  }

  for (i in  2:Sites){
    NatafdflistAutocorr[[i]]=NatafdflistAutocorr[[1]]
    SMAparams$GAFv[[i]]=SMAparams$GAFv[[1]]
    SMAparams$params[[i]]=SMAparams$params[[1]]
  }

  # Estimates parameters ai of SMA model for each site with FFT
  for (i in 1:1)
  {
    # Current Nataf transformeed GAF autoCOrr
    RGaf=SMAparams$GAFv[[i]]

    FFTLag=length(RGaf)-1
    g0=RGaf[1]
    g1=as.matrix(RGaf[-1])
    g2=flipdim(g1)

    co=c(g2,g0,g1)
    nM=2*FFTLag+1

    co_hat=(fft(z=co))

    ft=Re(ifft(x = sqrt(abs(co_hat))))
    ft1=ft[1]
    ftright=ft[2:(FFTLag+1)]
    ftleft=flipdim(matrix(ftright))

    ft0=c(ftleft,ft1,ftright)
    ft0=t(t(ft0))

    SMAparams$FFTa[[i]]=ft0
  }

  for (i in 2:Sites) {
    SMAparams$FFTa[[i]]= SMAparams$FFTa[[1]]
  }
  # correlation matrix of annual data between sites
  SMAparams$CorMat0=Cmat
  CCFupper=Cmat[upper.tri(x = Cmat, diag = F)]

  if (NatafIntMethod=='GH') {
    Nataf=NatafInvD(targetrho = CCFupper, fx = dist[i], fy = dist[i], paramlistfx = params[[i]], paramlistfy = params[[i]],
                    NatafIntMethod = 'GH', polydeg=polydeg, NoEval = NoEval, ...)
  } else if (NatafIntMethod=='Int') {
    Nataf=NatafInvD(targetrho = CCFupper, fx = dist[i], fy = dist[i], paramlistfx = params[[i]], paramlistfy = params[[i]],
                    NatafIntMethod = 'Int', polydeg=polydeg, NoEval = NoEval)
  } else if (NatafIntMethod=='MC') {
    Nataf=NatafInvD(targetrho = CCFupper, fx = dist[i], fy = dist[i], paramlistfx = params[[i]], paramlistfy = params[[i]],
                    NatafIntMethod = 'MC', polydeg=polydeg, NoEval = NoEval, ...)
  } else {
    print('Select an appropriate method to resolve Nataf integral. That is: "GH, "MC" or "Int"')
  }

  CorzMatrix=matrix(0, nrow=Sites, ncol=Sites)
  CorzMatrix[upper.tri(CorzMatrix)]=Nataf$rzEq
  CorzMatrix[lower.tri(CorzMatrix)]=rev(Nataf$rzEq)
  diag(CorzMatrix)=1
  SMAparams$CorMat=CorzMatrix

  # Perform NATAF transformation
  # We use this kind of loop to save time. Only the upper triangle is evaluated.
  # if (Sites!=1) {
  #   CorzMatrix=matrix(0,nrow=Sites,ncol=Sites)
  #   row.names(CorzMatrix)=namesSites
  #   colnames(CorzMatrix)=namesSites
  #   CCmat=SMAparams$CorMat0
  #   NatafdflistCrosscorr=list()
  #   CrossNames=NULL
  #   count=1
  #   for (i in 1:(Sites-1)) {
  #     for (j in (i+1):Sites){
  #
  #       paramslist$m1=params[[i]]
  #       paramslist$m2=params[[j]]
  #       rtarget=0
  #       icdf=c(dist[i],dist[j])
  #
  #       rtarget=CCmat[i,j]
  #
  #       if (NatafIntMethod=='GH') {
  #         Nataf=NatafInvD(targetrho = rtarget, fx = dist[i], fy = dist[j], paramlistfx = params[[i]], paramlistfy = params[[j]],
  #                         NatafIntMethod = 'GH', polydeg=polydeg, NoEval = NoEval, ...)
  #       } else if (NatafIntMethod=='Int') {
  #         Nataf=NatafInvD(targetrho = rtarget, fx = dist[i], fy = dist[j], paramlistfx = params[[i]], paramlistfy = params[[j]],
  #                         NatafIntMethod = 'Int', polydeg=polydeg, NoEval = NoEval)
  #       } else if (NatafIntMethod=='MC') {
  #         Nataf=NatafInvD(targetrho = rtarget, fx = dist[i], fy = dist[j], paramlistfx = params[[i]], paramlistfy = params[[j]],
  #                         NatafIntMethod = 'MC', polydeg=polydeg, NoEval = NoEval, ...)
  #       } else {
  #         print('Select an appropriate method to resolve Nataf integral. That is: "GH, "MC" or "Int"')
  #       }
  #       NatafdflistCrosscorr[[count]]=Nataf$dfnataf
  #       CrossNames=c(CrossNames,paste0('Site_',i,'Site_',j))
  #       CorzMatrix[i,j]=Nataf$rzEq
  #       count=count+1
  #     }
  #   }
  #   CorzMatrix=t(CorzMatrix)+CorzMatrix
  #   diag(CorzMatrix)=1
  #   SMAparams$CorMat=CorzMatrix
  # }


  # Find matrix of sums of ai parameters (used in the denominator of C)
  SSumAmatrix = matrix(0,nrow=Sites,ncol=Sites)
  row.names(SSumAmatrix)=namesSites
  colnames(SSumAmatrix)=namesSites
  for (i in 1:Sites)
  {
    for (j in 1:Sites)
    {
      SSumAmatrix[i,j]=sum(SMAparams$FFTa[[i]]*SMAparams$FFTa[[j]])
    }
  }
  SMAparams$SSumA=SSumAmatrix

  # Calcluate matrix "C"
  SMAparams$C=SMAparams$CorMat/SSumAmatrix
  SMAparams$C=round(SMAparams$C,3)
  # print(param$C)
  if (Sites!=1)
  { # for more than 1 sites, decomposition of the C matrix is needed
    # Decomposition of "C" matrix to derive parameter b(t)
    if(is.positive.definite(SMAparams$C))
    { # cholensky for positive definite matrices
      SMAparams$b=t(chol(SMAparams$C))
    } else
    { # smoothing of non-positive definite matrices to obtain positive
      if(DecoMethod=="cor.smooth") # method from "psych" package
      {
        tempC=cor.smooth(SMAparams$C)
      } else if (DecoMethod=="nearPD") # method from "Matrix" package
      {
        tempC=nearPD(SMAparams$C,corr=FALSE,keepDiag=TRUE)
        tempC=as.matrix(tempC$mat)
      }
      SMAparams$b=t(chol(tempC))
    }
  } else {
    SMAparams$b=1
  }
  names(NatafdflistAutocorr)=paste0('Site_',1:Sites)
  if (Sites>1) {
  # names(NatafdflistCrosscorr)=CrossNames
  # SMAparams$dfNatafCross=NatafdflistCrosscorr
  }
  SMAparams$dfNatafAuto=NatafdflistAutocorr


  return(SMAparams)
}
