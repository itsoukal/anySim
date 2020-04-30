#' @title Disaggregation of a coarse level timeseries sequence to a finer level timeseries sequence exhibiting the target marginal distribution and correlation structure (stationary).
#'
#' @description Disaggregation of a coarse level timeseries sequence to a finer level timeseries sequence exhibiting the target marginal distribution and correlation structure (stationary).
#'
#' @param HLSeries A vector, thati is coarse level timeseries sequence specifying the values to disaggreate into a time series sequence of finer level..
#' @param ARTApar A list containing the parameters of the model. The list is constructed by the function "EstARTAp".
#' @param max.iter A scalar specifying the maximum number of allowed repetitions (parameter of the disaggregation algorithm - typically set between 300-500.).
#' @param steps A scalar specifying the number of timesteps of the sequence to generate.
#' @param Adjust A logical operator (TRUE or FALSE) specifying whether (TRUE) or not (FALSE) to perfom the proportianal adjusting operation  (parameter of the disaggregation algorithm - typically set to TRUE).
#'
#' @return A list of the 3 generated time series (in vector format):
#' X: The final time series at the actual domain with the target marginal distribution and
#' correlation structure;
#'
#' @export
#'
#'
#' @examples
#' ## Disaggregation of a sequence of 500 steps of daily rainfall to 10-minute amounts
#' ## The lower-level process (i.e., that of 10-minute) is a zero-inflated one with p0=0.95, 
#' ## and a Burr Type XII distribution for the continuous part.
#' ## The target autocorrelation strucure is from CAS ACS.
#' \dontrun{
#' 
#' set.seed(124)
#' 
#' # Define the target autocorrelation structure of the finer-level process.
#' ACS=csCAS(param=c(1.688,1),lag=24) # CAS with b=1.688 and k=1.
#' 
#' # Define the target distribution function (ICDF).
#' FX='qmixed'# Define that distribution is of zero-inflated type.
#' 
#' # Define the distribution for the continuous part of the process.
#' # Here, a re-parameterized version of Burr Type-XII distribution is used.
#' qburr=function(p,scale,shape1,shape2) {
#'       require(ExtDist)
#'       x=ExtDist::qBurr(p=p,b=scale,g=shape1,s=shape2)
#'       return(x)
#' }
#' 
#' # Define the parameters of the zero-inflated distribution function.
#' pFX=list(p0=0.96,Distr=qburr,scale=0.181,shape1=7.642,shape2=0.296)
#'  
#' # Estimate the parameters of the auxiliary Gaussian AR(p) model.
#' param=EstARTAp(ACF=ACS,dist=FX,params=pFX,NatafIntMethod='GH')
#'
#' # Compose the daily series to be disaggregated
#' Sim=SimARTAp(ARTApar=param,burn=1000,steps=(24*6*500)) # generation of 10-min series
#' DailySeries=apply(X=matrix(data=Sim$X,ncol=24*6,byrow=1),MARGIN=1,FUN=sum)
#' 
#' ## Disaggregate daily series to 10-min data
#' disag10min=Disagg_ARTAp(HLSeries=DailySeries,ARTApar=param,max.iter=500,steps=24*6,Adjust=1)
#' }
Disagg_ARTAp<-function(HLSeries, ARTApar, max.iter, steps, Adjust=T){

  Zall=Zprevious=rnorm(length(ARTApar$phi))
  disag=matrix(NA, ncol=steps, nrow=length(HLSeries))

  for (i in 1:length(HLSeries)) {
    HLValue=as.vector(HLSeries[i])
    temp=Disagg_help(HLValue = HLValue, Zprevious = Zprevious, ARTApar = ARTApar,
                     max.iter = max.iter, steps = steps)
    disag[i,]=temp$X
    Zall=c(Zall, temp$Z)
    Zprevious=Zall[-c(1:(i*steps))]
    print(i)

  }
  X=as.vector(t(disag))
  return(list('X'=X))
}
