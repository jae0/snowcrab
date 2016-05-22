  map.set.information.diff = function(p, outdir, method="levelplot" ) {
    
    set = snowcrab.db( DS="set.complete")
    #variables = variable.list.expand("all.data")
    #variables = c('totmass.male.com', 'totmass.female.mat')
    variables = 'R0.mass'
    
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
      browser()
      
      # define compact list of variable year combinations for parallel processing
      mapyears = sort( unique(set$yr) )
      po = make.list( list(variables, mapyears ), Y=p )

      #for (i in p$init.files ) source( i )
     
        id <- c(1:po$nruns)
        id = as.numeric(id)

      for (i in id ) {
        v = po$runs[i,1]
        y = po$runs[i,2]
        outfn = paste( "test", sep=".")
        outloc = file.path( getwd() )
        xyz = set[ which(set$yr==y), c("plon","plat",v) ]
        xyz2 = set[ which(set$yr==y-1), c("plon","plat",v) ]
        diff <-xyz - xyz2
        er = empirical.ranges( db="snowcrab", v )  # range of all years
        datarange = seq( er[1], er[2], length.out=50)
        corners = data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))
        cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
        # cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
        names( diff) = c("plon", "plat", "z")
        
        map( diff, xyz.coords="planar", cfa.regions=T, depthcontours=T, pts=xyz, annot=y, fn=outfn, loc=outloc, at=datarange , col.regions=cols(length(datarange)+1), colpts=T, corners=p$planar.corners )
      }

    }
          
    return("Done")
  }


