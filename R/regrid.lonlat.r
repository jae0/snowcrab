
 
  regrid.lonlat = function (old, res, vr.to.sum=NULL, vr.to.cp=NULL) {
    old$gridid = paste(old$plon%/%res*res, old$plat%/%res*res, sep="." )
    new=data.frame( gridid = I(sort(unique(old$gridid)) ))
    for (vr in vr.to.sum) {
      l = aggregate(old[,vr], by=list(old$gridid), FUN=sum, na.rm=T)
      names(l) = c("gridid", vr)
      l$gridid = as.character(l$gridid)
      new = merge(x=new, y=l, by="gridid", all.x=T, all.y=F, sort=F)
    }
    if (!is.null(vr.to.cp) ) {
    for (vr in vr.to.cp) {
      l = aggregate(old[,vr], by=list(old$gridid), FUN=function(x) sort(unique(x))[1] )
      names(l) = c("gridid", vr)
      l$gridid = as.character(l$gridid)
      new = merge(x=new, y=l, by="gridid", all.x=T, all.y=F, sort=F)
    }}
    new$gridid = as.character( new$gridid )
    
    return( new )
  }



