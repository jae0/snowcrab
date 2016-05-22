
compute.proportions = function( x, var="egg" ) {
  
  x$dummy = 1
  ncrab=tapply( x$dummy, x$timevalue, FUN=length )    # number of primi crab on each day
  ncrab = as.data.frame( ncrab )
  ncrab$timevalue = as.numeric( rownames( ncrab))
  sub1=which( x[,var]==1 )              
  ncrab1=tapply( x$dummy[sub1], x$timevalue[sub1], FUN=length )
  ncrab1= as.data.frame (ncrab1)
  ncrab1$timevalue = as.numeric (rownames(ncrab1))
  prop = merge( ncrab, ncrab1, by=c("timevalue"), all.x=T, all.y=F) # combine total number of crab with var==1
  prop[is.na(prop)]=0   # turn na's to zero
  prop$fraction=NA
  prop$fraction= prop$ncrab1/prop$ncrab 
  return ( prop)

} 


