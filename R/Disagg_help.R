#' @title Disaggregation of a single value to a timeseries sequence exhibiting the target marginal distribution and correlation structure (stationary).
#'
#' @description Disaggregation of a single value to a timeseries sequence exhibiting the target marginal distribution and correlation structure (stationary).
#'
#' @param HLValue A scalar specifying the (single) value to disaggreate into a time series sequence.
#' @param ARTApar A list containing the parameters of the model. The list is constructed by the function "EstARTAp".
#' @param max.iter A scalar specifying the maximum number of allowed repetitions (parameter of the disaggregation algorithm - typically set between 300-500.).
#' @param steps A scalar specifying the number of timesteps of the sequence to generate.
#' @param Adjust A logical operator (TRUE or FALSE) specifying whether (TRUE) or not (FALSE) to perfom the proportianal adjusting operation  (parameter of the disaggregation algorithm - typically set to TRUE).
#'
#' @return A list of the 3 generated time series (in vector format):
#' X: The final time series at the actual domain with the target marginal distribution and correlation structure;
#' Z: The auxiliary Gaussian time series at the Gaussian domain and;
#' Diff: The difference between the (single) value to disaggreate with the sum of the generated sequence (i.e., aggregated to the same temporal level with the target value).
#' Dicrit: The sum of the generated sequence (i.e., aggregated to the same temporal level with the target value).
#'
#' @export
#'
#'
#' @examples
#' ## Disaggregation of a single 24-hour Rainfall ammount to 1-minute sequence.
#' ## The lower level process (i.e., that of 1-minute) is assumed to be a zero-inflated one
#' ## (i.e., using a mixed marginal distribution), with p0=0.95, and a Gamma distribution
#' ## for the continuous part with shape=1 and scale=8.
#' ## In this case, the target autocorrelation strucure is from
#' ## the CAS ACS with b=0 and k=0.7.
#' \dontrun{
#' FX='qmixed'
#' PFX=list(p0=0.95, Distr=qgamma, shape=1, scale=.8)
#' maxlag=60
#' ACFT=acsCAS(param = c(0, 0.2), lag=maxlag)
#' param=EstARTAp(ACF = ACFT, dist = FX, params = PFX, NatafIntMethod = 'GH')
#'
#' Sim=SimARTAp(ARTApar = param, burn = 1000, steps = 24*60*100) #24h*60min*100Days
#'
#' SUMX=apply(X = matrix(data = Sim$X, ncol=24*60, byrow = 0), MARGIN = 1, FUN = sum)
#' HLValue=as.vector(quantile(SUMX, probs = 0.5, type = 0));HLValue
# '
#' disag=Disagg_SV_S(HLValue = HLValue, ARTApar = param,
#'                           max.iter = 300, steps = 24*60, Adjust = 1)$X
#' }
Disagg_help<-function(HLValue, Zprevious, ARTApar, max.iter, steps, Adjust=T){
  
  phi=ARTApar$phi
  sigma=ARTApar$simga
  PFX=ARTApar$params
  FX=ARTApar$dist
  Zlist<-vector("list",max.iter)
  Xlist<-vector("list",max.iter)
  
  Zall<-rep(NA,steps)
  Distance<-rep(NA,max.iter)
  XBest=rep(NaN, steps)
  
  Dicrit=20000
  Thres=0.00002
  if (HLValue>Thres) {
    # while  ( data.table::between(Dicrit, HLValue*0.80, HLValue*1.20)==FALSE ) {
    while ( any(is.nan(XBest))==TRUE ) {
      Zall<-rep(NA,steps)
      for (j in 1:max.iter) {
        
        Zprevioustemp<-Zprevious
        Znew<-sum(rev(Zprevioustemp)*phi)+rnorm(1,mean=0,sd=sigma)
        Zall[1]<-Znew
        
        Zprevioustemp<-c(Zprevioustemp,Znew)
        
        for (index in 2:steps){ # need
          
          Zprevioustemp<-Zprevioustemp[-1]
          Znew<-sum(rev(Zprevioustemp)*phi)+rnorm(1,mean=0,sd=sigma)
          Zprevioustemp<-c(Zprevioustemp,Znew)
          Zall[index]<-Znew
          
        }
        U=pnorm(Zall, mean =0, sd = 1)
        Xall<-eval(as.call(c(as.name(FX),list(U),PFX)))
        XSum=sum(Xall,na.rm=TRUE)
        
        if (XSum==0) {
          Distance[j]=NA
        } else {
          Distance[j]=((HLValue-XSum))^2
        }
        
        Zlist[[j]]<-Zall
        Xlist[[j]]<-Xall
        
      }
      
      MinIndex=which.min(Distance)[1]
      
      # print(Xlist[[MinIndex]])
      ZBest=Zlist[[MinIndex]]
      XBest=Xlist[[MinIndex]]
      # Dicrit=abs(HLValue-sum(XBest))
      Dicrit=sum(XBest)
      
      Diff=HLValue/sum(XBest)
      
      if (is.infinite(Diff)==1) {
        XBest=rep(NaN,steps)
        max.iter=max(max.iter*5, 15000)
      }
      
      if (Adjust==T){ XBest=XBest*Diff }
      
    } # end while
  } else if (HLValue>0 & HLValue<=Thres){
    Diff=NA
    Dicrit=HLValue
    INDEX=sample(x = seq(1,steps), size = 1)
    XBest=rep(0,steps)
    XBest[INDEX]=HLValue
    # ZBest=rnorm(steps)
    Zprevioustemp<-Zprevious
    Znew<-sum(rev(Zprevioustemp)*phi)+rnorm(1,mean=0,sd=sigma)
    Zall[1]<-Znew
    Zprevioustemp<-c(Zprevioustemp,Znew)
    for (index in 2:steps){ # need
      Zprevioustemp<-Zprevioustemp[-1]
      Znew<-sum(rev(Zprevioustemp)*phi)+rnorm(1,mean=0,sd=sigma)
      Zprevioustemp<-c(Zprevioustemp,Znew)
      Zall[index]<-Znew
    }
    ZBest=Zall
    
  } else if (HLValue==0) {
    Diff=NA
    Dicrit=0
    XBest=rep(0,steps)
    # ZBest=rnorm(steps)
    Zprevioustemp<-Zprevious
    Znew<-sum(rev(Zprevioustemp)*phi)+rnorm(1,mean=0,sd=sigma)
    Zall[1]<-Znew
    Zprevioustemp<-c(Zprevioustemp,Znew)
    for (index in 2:steps){ # need
      Zprevioustemp<-Zprevioustemp[-1]
      Znew<-sum(rev(Zprevioustemp)*phi)+rnorm(1,mean=0,sd=sigma)
      Zprevioustemp<-c(Zprevioustemp,Znew)
      Zall[index]<-Znew
    }
    ZBest=Zall
  }
  
  return(list('X'=XBest, 'Z'=ZBest,'Diff'=Diff, 'Dicrit'=Dicrit))
}
