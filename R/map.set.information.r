  map.set.information = function(p, outdir, method="gmt" ) {
    
    set = snowcrab.db( DS="set.complete")
    variables = variable.list.expand("all.data")
    #variables = c('totmass.male.com', 'totmass.female.mat')
    
    if (method =="gmt") {
      # overrides to defaults
      p$tension = "-T.4"  # 0.35+ for steep; 0.25 for smooth
      p$maskres = "-S16k"
      p$interpres = "-nb"
      outdir = file.path( outdir, p$spatial.domain )
      gmt.map.variables( set, p=p, variables=variables, plottimes=p$plottimes, 
        basedir=outdir, conversions=p$conversions, init.files=p$init.files )
    }
    
    if (method=="levelplot") {
      
      # define compact list of variable year combinations for parallel processing
      mapyears = sort( unique(set$yr) )
      p = make.list( list(variables, mapyears ), Y=p )

      for (i in p$init.files ) source( i )
     
      if ( is.null(id)) id = c(1: p$nruns ) 
      id = as.numeric(id)

      for (i in id ) {
        v = p$runs[i,1]
        y = p$runs[i,2]
        outfn = paste( "test", sep=".")
        outloc = file.path( getwd() )
        xyz = set[ which(set$yr==y), c("plon","plat",v) ]
        er = empirical.ranges( db="snowcrab", v )  # range of all years
        datarange = seq( er[1], er[2], length.out=50)
        corners = data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))
        cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
        # cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
        names( xyz) = c("plon", "plat", "z")
        
        map( xyz, xyz.coords="planar", cfa.regions=T, depthcontours=T, pts=xyz, annot=y, fn=outfn, loc=outloc, at=datarange , col.regions=cols(length(datarange)+1), colpts=T, corners=p$planar.corners )
      }

    }
          
    return("Done")
  }


