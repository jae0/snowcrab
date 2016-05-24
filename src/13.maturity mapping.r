#maturity mapping

#loadfunctions('snowcrab','functionname="initialise.local.environment.r')
#loadfunctions('habitat')

di = file.path(project.datadirectory('snowcrab'),'R','maturity')
#gam model
load(file.path(di,'mat.gam.rdata'))
summary(.m4.gam)

 load(file.path(di,'x.mat.2014.rdata'))
 summary(x)

 p= list()

p$init.files = loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
p$libs = RLibrary( "mgcv", "chron", "lattice", "lattice", "grid", "fields", "parallel"  ) 

  
  p$regions = c("cfa4x", "cfanorth","cfasouth" )
  p$years.to.model=2004:2014

      p$habitat.threshold.quantile = 0.05 # quantile at which to consider zero-valued abundance
    	p$prediction.dyear = 9/12 # predict for ~ Sept 1 
      p$threshold.distance = 15  # limit to extrapolation/interpolation in km
     
 require(dplyr)   

params=p
   
hab = NULL
for(y in p$years.to.model) {
   PS = habitat.db ( DS="complete", year=y, p=params )
   PS$dyear = p$prediction.dyear
   PS$t = NA
         
        PST = temperature.interpolations( p=params, DS="spatial.interpolation", yr=y  )
				if (is.null(PST)) next ()
				  
        dyears = (c(1:(p$nw+1))-1)  / p$nw # intervals of decimal years... fractional year breaks
        dyr = as.numeric( cut( p$prediction.dyear, breaks=dyears, include.lowest=T, ordered_result=TRUE ) ) # integer representation of season
     
        PS$t = PST[, dyr ]
        PS$t[ which(PS$t < -2) ] = -2
			  PS$t[ which(PS$t > 30) ] = 30 
	    PS = PS[,c('yr','z','t','plon','plat')]
	    PS$grd = 1:nrow(PS)
	    cw = seq(40,126,by=2)
	   	 io = dim(PS)
	    PPS = cbind(PS[rep(seq_len(io[1]),each=length(cw)),],cw=rep(cw,times=io[1]))
        re = predict(.m4.gam,newdata=PPS,type='response')
        attr(re,'dimnames') <- NULL
        re = c(re)
        PPS = cbind(PPS,re=re)
      	save(PPS,file=file.path(di,paste(y,'complete.data.mat50.rdata')))
       }     

#summarize for plot
for(y in p$years.to.model) {
	load(file=file.path(di,paste(y,'complete.data.mat50.rdata')))
      	pp = group_by(PPS,grd)
      	ap = summarize(pp,ld50=findMax(cw,abs(re-0.5)*-1)[1])
      out = unique(pp[,c('yr','plon','plat','grd')])
      out = merge(out,ap,by='grd')
        outfn = paste( y,"test", sep=".")
        outloc = file.path( getwd() )
        xyz = out[, c("plon","plat",'ld50') ]
        xyz$ld50[which(xyz$ld50<50 | xyz$ld50>150)]<- NA
        er = range(xyz$ld50,na.rm=T)
        datarange = seq( er[1], er[2], length.out=50)
        corners = data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))
        cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
        # cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
        names( xyz) = c("plon", "plat", "z")
        
        map( xyz, xyz.coords="planar", cfa.regions=T, depthcontours=T,  annot=y, fn=outfn, loc=outloc, at=datarange , col.regions=cols(length(datarange)+1), colpts=T, corners=corners )
}
