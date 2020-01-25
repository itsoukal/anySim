Disagg_Seq_S<-function(HLSeries, ARTApar, max.iter, steps, Adjust=T){
  
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