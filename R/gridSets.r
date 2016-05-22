gridSets <- function(dat,grid.size=50) {
 
  print('grid size is set for km only')
 
  dat <- dat[complete.cases(dat[,c('lon','lat')]),]

  xwin <- data.frame(lat=c(min(dat$lat),min(dat$lat),max(dat$lat),max(dat$lat)),lon=c(min(dat$lon),max(dat$lon),min(dat$lon),max(dat$lon)))
 
  lon.d <- vincenty(loc1=xwin[1,],loc2=xwin[2,],a=6378.13700,f=1/298.257223563)
  lat.d <- vincenty(loc1=xwin[1,],loc2=xwin[3,],a=6378.13700,f=1/298.257223563)
 
  dplo <- (xwin[2,2]-xwin[1,2])/lon.d # deg NS per km
  dpla <- (xwin[3,1]-xwin[2,1])/lat.d # deg EW per km
  
  seqlo <- dplo*grid.size
  seqla <- dpla*grid.size             
   
 
 ys <- seq(xwin[1,1],xwin[3,1],by=seqla)
 xs <- seq(xwin[1,2],xwin[2,2],by=seqlo)
 require(PBSmapping)
 gr <- makeGrid(x=xs,y=ys,projection="LL",zone=20,addSID=F)
 
 #make pbsmap friendly
	dat$X 	<- dat$lon
	dat$Y 	<- dat$lat
	dat$EID <- 1:nrow(dat)
	attr(dat,'projection') <- "LL"	
	fp 		<- findPolys(dat,gr,maxRows=1e6)	
 	dat <- merge(dat,fp,by='EID')
  return(list(dat,gr))
}