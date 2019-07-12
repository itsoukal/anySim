#' @title Estimation of the auxiliary PAR(1) model parameters
#'
#' @description Estimation of parameters of the auxiliary PAR(1) model to simulate the auxiliary cyclostationary Gaussian process to establish the target season-to-season correlation structure.
#'
#' @param s2srtarget A k-dimensional vector with the lag-1 season-to-season correlation coefficients.
#' @param dist A k-dimensional string vector indicating the quantile function of the target marginal distributions (i.e., the ICDFs).
#' @param params A k-dimensional named list with the parameters of the target distributions.
#' @param NatafIntMethod A string ("GH", "Int", or "MC") indicating the integration method to resolve the Nataf integral.
#' @param NoEval A scalar indicating (default: 9) the number of evaluation points for the integration methods.
#' @param polydeg A scalar indicating the order of the fitted polynomial. If polydeg=0, then another curve is fitted.
#' @param ... Additional named arguments for the selected "NatafIntMethod" method.
#'
#' @note Avoid the use of the "GH" method (i.e., NatafIntMethod='GH'), when the marginal(s) are discrete.
#'
#' @return A list with the parameters of the cyclostationary auxiliary Gaussian PAR(1) model.
#' @export
#'
#' @examples
#' ## Parameter estimation of a cyclostationary process with 12 seasons.
#' \dontrun{
#' rtarget<-c(0.5, 0.7, 0.6, 0.4, 0.5, 0.7, 0.8, 0.7, 0.6, 0.4, 0.5, 0.7)
#' NumOfSeasons=length(rtarget)
#'
#' FXs<-rep('qmixed',NumOfSeasons)
#' PFXs<-vector("list",NumOfSeasons)
#' PFXs<-lapply(PFXs,function(x) x<-list(p0=0.4,Distr=qexp,rate=0.5))
#'
#' SPARTApar<-EstSPARTA(s2srtarget = rtarget, dist = FXs,
#' params = PFXs, NatafIntMethod = 'GH', NoEval = 9, polydeg = 8)
#' }

EstSPARTA=function(s2srtarget, dist, params, NatafIntMethod='GH', NoEval=9, polydeg=6, ...){

    NumOfSeasons<-length(dist)

    r<-rep(NA,NumOfSeasons)
    Natafdflist=list()

    ft<-dist[1]
    ftt<-dist[NumOfSeasons]
    pft<-params[[1]]
    pftt<-params[[NumOfSeasons]]

    # rz<-seq(0,1,length.out = NoEval)
    # rx<-NatafGH(rz,ft,ftt,pft,pftt)
    #
    # r[1]<-NatafInv(rz,rx,target = s2srtarget[1],polydeg=polydeg)


    if (NatafIntMethod=='GH') {
      Nataf=NatafInvD(targetrho = s2srtarget[1], fx = ft, fy = ftt, paramlistfx = pft, paramlistfy = pftt,
                      NatafIntMethod = 'GH', polydeg=polydeg, NoEval = NoEval, ...)
    } else if (NatafIntMethod=='Int') {
      Nataf=NatafInvD(targetrho = s2srtarget[1], fx = ft, fy = ftt, paramlistfx = pft, paramlistfy = pftt,
                      NatafIntMethod = 'Int', polydeg=polydeg, NoEval = NoEval)
    } else if (NatafIntMethod=='MC') {
      Nataf=NatafInvD(targetrho = s2srtarget[1], fx = ft, fy = ftt, paramlistfx = pft, paramlistfy = pftt,
                      NatafIntMethod = 'MC', polydeg=polydeg, NoEval = NoEval, ...)
    } else {
      print('Select an appropriate method to resolve Nataf integral. That is: "GH, "MC" or "Int"')
    }
    Natafdflist[[1]]=Nataf$dfnataf
    r[1]<-Nataf$rzEq

    for(i in 2:NumOfSeasons){

      ft<-dist[i]
      ftt<-dist[(i-1)]
      pft<-params[[i]]
      pftt<-params[[i-1]]

      # rz<-seq(0,1,length.out = NoEval)
      # rx<-NatafGH(rz,ft,ftt,pft,pftt)
      #
      # r[i]<-NatafInv(rz,rx,target = s2srtarget[i],polydeg=polydeg)
      if (NatafIntMethod=='GH') {
        Nataf=NatafInvD(targetrho = s2srtarget[i], fx = ft, fy = ftt, paramlistfx = pft, paramlistfy = pftt,
                        NatafIntMethod = 'GH', polydeg=polydeg, NoEval = NoEval, ...)
      } else if (NatafIntMethod=='Int') {
        Nataf=NatafInvD(targetrho = s2srtarget[i], fx = ft, fy = ftt, paramlistfx = pft, paramlistfy = pftt,
                        NatafIntMethod = 'Int', polydeg=polydeg, NoEval = NoEval)
      } else if (NatafIntMethod=='MC') {
        Nataf=NatafInvD(targetrho = s2srtarget[i], fx = ft, fy = ftt, paramlistfx = pft, paramlistfy = pftt,
                        NatafIntMethod = 'MC', polydeg=polydeg, NoEval = NoEval, ...)
      } else {
        print('Select an appropriate method to resolve Nataf integral. That is: "GH, "MC" or "Int"')
      }
      Natafdflist[[i]]=Nataf$dfnataf
      r[i]<-Nataf$rzEq


    }
    names(r)=paste0('Season_', 1:NumOfSeasons)
    names(Natafdflist)=paste0('Season_', 1:NumOfSeasons)
    names(dist)=paste0('Season_', 1:NumOfSeasons)
    names(params)=paste0('Season_', 1:NumOfSeasons)

    return(list('NumOfSeasons'=NumOfSeasons,'s2sEq'=r, 'Natafdf'=Natafdflist,'dist'=dist,'params'=params))
  }
