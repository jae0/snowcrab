map.fisheries.data = function(p, outdir,  FUN, yrs, variable='effort',probs=c(0,0.975),log.variable=F,offset,theta=75,pres=10) {

  x = logbook.db( DS="logbook" )

  x = x [polygon_inside( x, region="isobath1000m"),]
  x = x[ which(x$effort <= 300) ,]
  x = x[ which(x$cpue < 500),]
  x$year=x$yr #this creates proper offset for 4X, 2017-18 season =2017
  x = lonlat2planar( x,  proj.type=p$aegis_proj4string_planar_km )
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
    xd$plon = grid_internal( xd$plon, g$plon )
    xd$plat = grid_internal( xd$plat, g$plat )
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
    dir.create (outloc, showWarnings=FALSE, recursive =TRUE)
    fn = file.path( outloc, paste(outfn, "png", sep="." ) )
    png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
      lp = aegis_map( xyz[[i]], xyz.coords="planar", depthcontours=TRUE,
        pts=NULL, annot=yrs[i], annot.cex=4, at=datarange, colorkey=NULL,
        col.regions=cols(length(datarange)+1), colpts=FALSE, corners=p$corners, display=FALSE, rez=c(pres,pres), plotlines="cfa.regions")
    print(lp)
    dev.off()
    print(fn)


    #try( ) ##BZ Nov 2018- Removed. Producing errors. No expression to "try" within the try() function
  }
}
