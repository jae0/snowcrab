
  snowcrab.habitat.db = function( ip=NULL, DS="", p, v=NULL, y=NULL ) {

    outdir = file.path( project.datadirectory("snowcrab"), "R", "gam", "habitat" )
    dir.create(path=outdir, recursive=T, showWarnings=F)
    
    if (DS=="PS" ) {
      out = NULL
        fn.PS = file.path( outdir, paste( "PS", v, y, "rdata", sep="." ) )
        if ( ! (file.exists( fn.PS)) ) return(NULL)
        load(fn.PS)
      return (PS)
    }
    
    if (DS=="K" ) {
      out = NULL
        fn.K = file.path( outdir, paste( "K", v, y, "rdata", sep="." ) )
        if ( ! (file.exists( fn.K)) ) return(NULL)
        load(fn.K)
      return (K)
    }
     
    if (!is.null(p$init.files)) for( i in p$init.files ) source (i)
    if (is.null(ip)) ip = 1:p$nruns

    for ( iip in ip ) {

        y = p$runs[iip,"y"]
        v = p$runs[iip,"v"]
    
        v0 = v = p$runs[iip,"v"]
        if ( v0 =="R0.mass.environmentals.only" ) v="R0.mass"


        PS = habitat.db ( DS="complete", year=y, p=p )
				PS$dyear = p$prediction.dyear  # must be same as above
				PS$t = NA
         
        PST = temperature.db( p=p, DS="spatial.interpolation", yr=y  )
				if (is.null(PST)) next ()
				    
        dyears = (c(1:(p$nw+1))-1)  / p$nw # intervals of decimal years... fractional year breaks
        dyr = as.numeric( cut( p$prediction.dyear, breaks=dyears, include.lowest=T, ordered_result=TRUE ) ) # integer representation of season
        PS$t = PST[, dyr ]
        PS$t[ which(PS$t < -2) ] = -2
			  PS$t[ which(PS$t > 30) ] = 30 

        iitna = which( ! is.finite( PS$t ) ) 
        if (length(iitna)>0) PS$t[iitna] = PS$tmean[iitna]

        PS$z = log(PS$z)
        PS$dt.seasonal = PS$tmean - PS$t 
        PS$dt.annual = PS$tmean - PS$tmean.cl
        PS$sa = 1
        
        if ( y < 1998) PS$yr = floor(median( p$years.to.model ))  # assume similar conditions as those found in 1998 for the year-effect (no extrapolation)

				# predictions
        # Alternate models using only environmental information without years
        Habitat.model =  habitat.model.db( DS="habitat", v=v0 )
        preds = predict( Habitat.model, newdata=PS, type="response", na.action="na.pass" , se.fit=TRUE, block.size=10000, newdata.guaranteed=T )
        
        PS$habitat.mean = preds$fit
        PS$habitat.sd = preds$se.fit
        rm (preds); gc()
        
        # levelplot( habitat.mean ~ plon + plat, data=PS, aspect="iso" )

        iHabitat = which( PS$habitat.mean > p$habitat.threshold.quantile   & (PS$habitat.mean - 2 * PS$habitat.sd) > 0  )

        # levelplot( habitat.mean ~ plon + plat, data=PS[iHabitat,], aspect="iso" )

				totalsurfacearea = length( iHabitat ) * (p$pres*p$pres)
        temp.mean = mean( PS$t, na.rm=T )
        temp.sd   = sd(   PS$t, na.rm=T )
        
        fn.res = file.path( outdir, paste( "PS", v0, y, "rdata", sep="." ) )
        print (fn.res )
        save( PS, file=fn.res, compress=T )

        # MAP

          datarange = seq( 0, 1, length.out=150)
          cols = color.code( "seis", datarange )
          outfn = paste( "prediction.habitat.mean.direct", v0, y, sep=".")
          map( xyz=PS[,c("plon", "plat", "habitat.mean")], cfa.regions=T, depthcontours=T, pts=NULL, annot=paste( v, y ), 
            annot.cex=p$annot.cex, corners=planar.corners, fn=outfn, loc=outdir, at=datarange, 
            col.regions=cols, rez=c(p$pres,p$pres) )

        K = NULL
        for (r in p$regions ){

          iRegion.PS = filter.region.polygon(x=PS[ , c("plon", "plat")], region=r, planar=T)
          iEastof250 = which( PS$plon > 250 )
          iRegion.PS = intersect( iRegion.PS, iEastof250 )
					iHabitatSubarea = intersect( iRegion.PS, iHabitat ) # plotting surface and real habitat area
          sa.region = length( iHabitatSubarea ) * (p$pres*p$pres)
              
          L = data.frame( yr=y, vars=v0, region=r, 
              sa.total=totalsurfacearea, sa.region=sa.region, 
              temp.total=temp.mean, temp.total.sd= temp.sd,
              temp.region= mean( PS$t[iHabitatSubarea], na.rm=T ), 
              temp.region.sd=sd( PS$t[iHabitatSubarea], na.rm=T ),
              datestamp=as.character( Sys.time() ) )
          K = rbind( K, L )    
        } # regions
        
        fn.K = file.path( outdir, paste( "K", v0, y, "rdata", sep="." ) )
        save( K, file=fn.K, compress=T )


    }
  
  }



