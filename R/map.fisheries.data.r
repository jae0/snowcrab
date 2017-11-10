map.fisheries.data = function(p, outdir,  FUN, yrs, variable='effort',probs=c(0,0.975),log.variable=F,offset,theta=75,pres=10) {

  x = logbook.db( DS="logbook" )
  
  x = x [emgis::polygon_inside( x, region="isobath1000m"),]
  x = x[ which(x$effort <= 300) ,]
  x = x[ which(x$cpue < 500),]
  
  x = lonlat2planar( x,  proj.type=p$internal.projection )
  x = subset(x,select=c('year','plon','plat',variable))
  names(x)[4] = 'z'
  
  if(variable=='landings')x$z = x$z/1000

  predlocs = bathymetry.db(p=p, DS="baseline")
  er = quantile( x$z[x$z>0], probs=probs)  # range of all years
  ler = er
  
  if(missing(offset))offset = er[1]  # offset fot log transformation
  if(missing(yrs))yrs=sort(unique(x$year))
        
  xyz = list()

  for(i in 1:length(yrs)){
  
    xd = subset(x,year==yrs[i],c('plon','plat','z'))
    
    p$pres=pres
    g= spatial_grid(p,"planar.coords")
    xd$plon = grid.internal( xd$plon, g$plon )
    xd$plat = grid.internal( xd$plat, g$plat )
    xyz[[i]] = aggregate(z~plon+plat,data=xd,FUN)
    names( xyz[[i]]) = c("plon", "plat", "z")
  }

  # calculate data range to show
  ler=quantile( do.call("rbind",xyz)$z, probs=probs)
  datarange = seq( ler[1], ler[2], length.out=50)
  #ckey=list(labels=list(at=pretty(datarange,8),labels=pretty(datarange,8),cex=2))

  cols = colorRampPalette(c( "yellow", "orange","darkred", "black"), space = "Lab")
  #cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
  for(i in 1:length(yrs)){
    xyz[[i]]$z[xyz[[i]]$z>ler[2]] = ler[2] # set levels higher than max datarange to max datarange
    outfn = paste( variable,yrs[i], sep=".")
    outloc = file.path( outdir,variable)

    # do the map
    try( emei::emei_map( xyz[[i]], xyz.coords="planar", cfa.regions=T, depthcontours=T, pts=NULL, annot=yrs[i], fn=outfn, loc=outloc, at=datarange,colorkey=NULL , col.regions=cols(length(datarange)+1), colpts=F, corners=p$corners, display=F,rez=c(pres,pres)))
  }
}


