


#TODO BC add functionality for pdf&kml outputs

map.set.information = function(p, outdir, variables, mapyears, interpolate.method='tps', theta=p$pres*25, ptheta=theta/2.3,
                               idp=2, log.variable=TRUE, add.zeros=TRUE, minN=10, probs=c(0.025, 0.975) ) {

    set = snowcrab.db( DS="set.biologicals")
    if(missing(variables)){
      variables = bio.snowcrab::snowcrab.variablelist("all.data")
      variables = intersect( variables, names(set) )
    }

    # define compact list of variable year combinations for parallel processing
    if (missing(mapyears)) mapyears = sort( unique(set$yr) )
    if (exists( "libs", p)) RLibrary( p$libs )

    predlocs = spatial_grid(p) 
    predlocs = planar2lonlat( predlocs,  proj.type=p$aegis_proj4string_planar_km   )
    predlocs$z = bathymetry_lookup_rawdata( spatial_domain=p$spatial_domain, M=predlocs[, c("lon", "lat")] ) 
    predlocs = predlocs[ geo_subset( spatial_domain=p$spatial_domain, Z=predlocs ), ]
  
    for ( v in variables ) {
      for ( y in mapyears ) {
     
          ratio=FALSE
          outfn = paste( v,y, sep=".")
          outloc = file.path( outdir,v)
          ref=y

          if (grepl('ratio', v)) ratio=TRUE

          set_xyz = set[ which(set$yr==y), c("plon","plat",v) ]
          names( set_xyz) = c("plon", "plat", "z")
          set_xyz = na.omit(subset(set_xyz,!duplicated(paste(plon,plat))))
          if(nrow(set_xyz)<minN)next() #skip to next variable if not enough data


          offset = empirical.ranges( db="snowcrab", v, remove.zeros=T , probs=0)  # offset fot log transformation
          er = empirical.ranges( db="snowcrab", v, remove.zeros=T , probs=probs)  # range of all years
          if(ratio)er=c(0,1)
          ler = er
          withdata=which(set_xyz$z > 0)
          if (length(withdata) < 3) print(paste("skipped",v, y, "<3 data points to create map", sep="."))
          if (length(withdata) < 3) next()
          S = set_xyz[ withdata, c("plon", "plat") ]


          distances =  rdist( predlocs[,c("plon", "plat")], S)
          distances[ which(distances < ptheta) ] = NA
          ips = which( !is.finite( rowSums(distances) ) )
          plocs=predlocs[ips,]
          if(log.variable){
            set_xyz$z = log(set_xyz$z+offset)
            ler=log(er+offset)
            #if(offset<1)if(shift) xyz$z = xyz$z + abs(log(offset))
          }

          datarange = seq( ler[1], ler[2], length.out=50)
  #
          #if(logit.variable){
          #  sr=set[,v]
          #  sr=sr[sr>0&sr<1&!is.na(sr)]
          #  lh=range(sr)
          #  set_xyz$z[set_xyz$z==0] = lh[1]
          #  set_xyz$z[set_xyz$z==1] = lh[2]
          #  set_xyz$z = logit(set_xyz$z)
          #  #if(offset<1)if(shift) xyz$z = xyz$z + abs(log(offset))
          #  ler=logit(quantile(sr,probs))
          #  datarange = seq( ler[1], ler[2], length.out=50)
          #}
          xyzi = na.omit(set_xyz)

          if(nrow(xyzi)<minN||is.na(er[1]))next() #skip to next variable if not enough data

          #!# because 0 in log space is actually 1 in real space, the next line adds the log of a small number (offset)
          #!# surrounding the data to mimic the effect of 0 beyond the range of the data
          if(add.zeros)  xyzi =na.omit( zeroInflate(set_xyz,corners=p$corners,type=2,type.scaler=0.5,eff=log(offset),blank.dist=20) )

          if(interpolate.method=='tps'){

            u= fastTps(x=xyzi[,c("plon","plat")] , Y=xyzi[,'z'], theta=theta )
            res = cbind( plocs[,1:2], predict(u, xnew=plocs[,1:2]))
          }
          if(interpolate.method=='idw'){
            require(gstat)
            u = gstat(id = "z", formula = z ~ 1, locations = ~ plon + plat, data = xyzi, set = list(idp = idp))
            res = predict(u, plocs[,1:2])[,1:3]
          }
          #print(summary(set_xyz))
          #print(summary(res))

          xyz = res
          names( xyz) = c("plon", "plat", "z")
          #if(shift)xyz$z = xyz$z - abs(log(offset))

          cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")

          xyz$z[xyz$z>ler[2]] = ler[2]
          if(ratio)xyz$z[xyz$z<ler[1]] = ler[1]

          ckey=NULL
          if(log.variable){
            # create labels for legend on the real scale
            labs=as.vector(c(1,2,5)%o%10^(-4:5))
            labs=labs[which(labs>er[1]&labs<er[2])]
            ckey=list(labels=list(at=log(labs+offset),labels=labs,cex=2))
          }

          dir.create (outloc, showWarnings=FALSE, recursive =TRUE)
          annot=ref
          filename=file.path(outloc, paste(outfn, "png", sep="."))
          print(filename)
          png( filename=filename, width=3072, height=2304, pointsize=40, res=300 )
          lp = aegis_map( xyz, xyz.coords="planar", depthcontours=TRUE, pts=set_xyz[,c("plon","plat")],
            annot=annot, annot.cex=4, at=datarange , col.regions=cols(length(datarange)+1),
            colpts=F, corners=p$corners, display=F, colorkey=ckey, plotlines="cfa.regions" )
          print(lp)
          dev.off()

      }
    }

    return("Done")
  }
