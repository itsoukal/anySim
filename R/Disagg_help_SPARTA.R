






Disagg_help_SPARTA=function(HLValue, Zprevious, SPARTApar, max.iter, steps, Adjust=T) {

  NumOfSeasons <- SPARTApar$NumOfSeasons
  FXs = SPARTApar$dist
  PFXs = SPARTApar$params
  r = SPARTApar$s2sEq

  Z = X = U = matrix(data = NA, nrow = max.iter, ncol = NumOfSeasons)
  m = nrow(Z)
  n = ncol(Z)

  rmod = (1 - r^2)^0.5

  for (i in 1:max.iter) {
    Zprevious_temp=Zprevious
    for (j in 1:n) {
      if (j != 1) {
        Z[i, j] = Z[i, j - 1] * r[j] + rmod[j] * rnorm(1)
      }
      else if (i == 1 && j == 1) {
        Z[i, j] = Zprevious_temp
      }
      else if (j == 1) {
        Z[i, j] = Zprevious_temp * r[j] + rmod[j] * rnorm(1)
      }
    }

    for (j in 1:NumOfSeasons) {
      U[i, j] = pnorm(Z[i, j])
      fs = FXs[j]
      pfs = PFXs[[j]]
      X[i, j] = eval(as.call(c(as.name(fs), list(U[i, j]), pfs)))
    }

  }

  Distance = log(abs(Matrix::rowSums(X)-rep(HLValue, m)))
  idbest=which.min(Distance)
  Diff=Distance[idbest]
  XBest=X[idbest,]
  ZBest=Z[idbest,]
  if (Adjust==T) {
    Diff_Eucl=HLValue-sum(XBest)
    XBest=X[idbest,]= XBest*HLValue/sum(XBest)
    ZBest=Z[idbest,]= ZBest
  }

  colnames(Z) = paste0("Season_", 1:NumOfSeasons)
  colnames(X) = paste0("Season_", 1:NumOfSeasons)
  colnames(U) = paste0("Season_", 1:NumOfSeasons)
  return(list('X'=XBest, 'Z'=ZBest,'Diff'=Diff,'Diff_Eucl'=Diff_Eucl))
}
