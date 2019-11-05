biomass.summary.table= function(x){

  y = res$mcmc
  sb= res$p$fishery_model$standata
  
  ntacs = sb$nProj
  yrs0 = as.numeric( as.character( rownames(sb$IOA) ) )
  yrs = c( yrs0, (max(yrs0)+c(1:sb$M) ) )
  
SI =  apply( y$q, 2, median, na.rm=T  )

tab=as.data.frame(yrs)
for (i in 1:3) {
  meanval = apply( y$B[,,i], 2, mean, na.rm=T  )
  meanval = as.data.frame(meanval)
  if (i==1) {names(meanval)="nens"}
  if (i==2) {names(meanval)="sens"}
  if (i==3) {names(meanval)="4x"}
  tab=cbind(tab, meanval) 
}
tab=tab[which(tab$yrs %in% yrs0),]
print("The following table is the mean modelled biomass estimates by area by year:")
print("Area names hard-coded, confirm they are logical")

print(tab)
}