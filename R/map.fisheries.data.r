
map.fisheries.data = function(p, outdir,  FUN, yrs, variable='effort',method='tps',probs=c(0,0.975),log.variable=T,offset,theta=75) {

  x = logbook.db( DS="logbook" )
  
  x = x [filter.region.polygon( x, region="isobath1000m"),]
  x = x[ which(x$effort <= 300) ,]
  x = x[ which(x$cpue < 500),]
  
  x = lonlat2planar( x,  proj.type=p$internal.projection )
  x = subset(x,select=c('year','plon','plat',variable))
  names(x)[4] = 'z'
  
  predlocs = bio.bathymetry::bathymetry.db(p=p, DS="baseline")
  er = quantile( x$z[x$z>0], probs=probs)  # range of all years
  ler = er
  
  if(missing(offset))offset = er[1]  # offset fot log transformation
  if(log.variable)ler=log(er+offset)
  if(missing(yrs))yrs=sort(unique(x$year))
        
  datarange = seq( ler[1], ler[2], length.out=50)
  
  for(i in 1:length(yrs)){
  
    outfn = paste( variable,yrs[i], sep=".")
    outloc = file.path( outdir,variable)
    xd = subset(x,year==yrs[i],c('plon','plat','z'))
    
    if(method=='tps'){
      browser()
      xd = aggregate(z~plon+plat,data=xd,FUN)
      if(log.variable)xd$z = log(xd$z+offset)
      u= fastTps(x=xd[,c("plon","plat")] , Y=xd[,'z'], theta=theta )
      res = cbind( predlocs[,1:2], predict(u, xnew=predlocs[,1:2]))
    }
    xyz = res
    names( xyz) = c("plon", "plat", "z")
    cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
    xyz$z[xyz$z>ler[2]] = ler[2]
    ckey=NULL

    if(log.variable){
      labs=as.vector(c(1,2,5)%o%10^(-4:5))
      labs=labs[which(labs>er[1]&labs<er[2])]
      ckey=list(labels=list(at=log(labs+offset),labels=labs,cex=2))
    }
    
    try( map( xyz, xyz.coords="planar", cfa.regions=T, depthcontours=T, pts=set_xyz[,c("plon","plat")], annot=yrs[i], fn=outfn, loc=outloc, at=datarange , col.regions=cols(length(datarange)+1), colpts=F, corners=p$corners, display=F,colorkey=ckey))
  }
}


