#interpolate catch rates from survey data for mapping

#####
#####Notes this was run to produce the interpolated survey maps for EMERA fixed link November 15, 2013 12:27:44 PM Maritime Link


p = bio.snowcrab::load.environment()
require(raster)

	  do.interpolation=F

if(do.interpolation) {
#get crab data new
	det <- snowcrab.db('set.biologicals')
	gps <- c('totno.male.imm','totno.female.imm','totno.male.mat','totno.female.mat')
	det <- det[,c('trip','station','set','plat','plon',gps,'yr')]
	d 	<- aggregate(trip~station,data=det,FUN=length)
	det	<- det[det$yr>2003 & det$station %in% d[d[,2]>3,1],]

	geo.mean=T
  
	for(i in 1:length(gps)) {
		de 				<- det[,c('station',gps[i],'plat','plon','trip','set')]
		names(de)[2] 	<- 'gg'

		if(!geo.mean) dp 	<- aggregate(cbind(gg,plat,plon)~station,data=de,FUN=mean)
		if(geo.mean)  dp	<- aggregate(cbind(gg,plat,plon)~station,data=de,FUN=function(x) exp(mean(log(x+1)))-1)

		dp$limm 		<- log(dp[,'gg']+1)
		params 			<- list(nmax=15,drange=15,power=0.8)
# 		fp 				<- interpol.grid(dp[,c('plon','plat','limm')],locs=spatial_grid(p=p, nm=c("x","y")) ,method='inv.dist.gstat',params=params,outfile='iminterp.dat')
    fp        <- iwd(dp[,c('plon','plat','limm')],locs=spatial_grid(p=p, nm=c("x","y")) ,method='inv.dist.gstat',params=params,outfile='iminterp.dat')
		er 				<- with(dp[dp$limm>0,],quantile(limm,probs=c(0.025,0.975)))
		datarange 		<- c(0,seq( er[1], er[2], length.out=100))
    corners 		<- data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))
    cols 			<- colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
    names(fp)[1:3] 	<- c('plon','plat','z')
    fo 				<- file.path("C:","~","bio","snowcrab","Adam","Emera","abundance maps","global average")
    if(geo.mean) fo 				<- file.path("C:","~","bio","snowcrab","Adam","Emera","abundance maps","global average","geomean")
    if(!geo.mean) fo 				<- file.path("C:","~","bio","snowcrab","Adam","Emera","abundance maps","global average","arithmean")

      outfn = gps[i]
      dir.create (fo, showWarnings=FALSE, recursive =TRUE)
      fn = file.path( fo, paste(outfn, "png", sep="." ) )
      png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
      lp = aegis::aegis_map( fp[,1:3], xyz.coords="planar", cfa.regions=T, depthcontours=T, 
        at=datarange , col.regions=cols(length(datarange)+1), colpts=T, corners=p$planar.corners, annot=gps[i] )
      print(lp)
      dev.off()

    ras 			<- rasterFromXYZ(fp[,c('plon','plat','z')],crs=p$internal.crs)

    writeRaster(ras,filename=file.path(fo,paste(gps[i],".asc",sep="")), format="ascii", overwrite=TRUE)
    }

}
