
  interpolation.db = function( ip=NULL, DS=NULL, p=NULL, 
    varnames = c("snowcrab.large.males_abundance", "snowcrab.large.males_presence_absence") ) {

    if (DS %in% c( "biomass", "biomass.redo" )) {
    
      outdir = file.path( project.datadirectory("bio.snowcrab"), "modelled", "biomass" )
      dir.create(path=outdir, recursive=T, showWarnings=F)
      fn = file.path( outdir, "biomass.rdata" )

      if (DS=="biomass") {
        B = NULL
        if ( file.exists( fn) ) load(fn)
        return (B)
      }

      
      set = bio.indicators::survey.db( p=p, DS="set.filter" ) # mature male > 95 mm 
 
      ii = which( set$totmass > 0 )
      qs = quantile( set$totmass[ii], probs=p$lbm_quantile_bounds, na.rm=TRUE )
      qs = log(qs) # 

      bm = snowcrab_lbm( p=p, DS="baseline", ret="mean", varnames=varnames )
      m = bm[[1]]  # biomass
      h = bm[[2]]  # habitat
      bm= NULL
      
      ll = which(h < p$habitat.threshold.quantile )
      if (length(ll) > 0 ) h[ll] = 0
      m = log( exp(m) * h ) # bm[[2]] is serving as weight/probabilities
      # #respect the bounds of input data (no extrapolation)
      
      qq = which( m < qs[1] )
      if (length(qq) > 0 ) m[qq] = NA
      rr = which( m > qs[2] )
      if (length(rr) > 0 ) m[rr] = qs[2]

      if(0) {
        bs = bio.bathymetry::bathymetry.db(p=p, DS="baseline")
        levelplot( m[,16] ~ plon+plat, bs, aspect="iso")
        for (i in 1:16) print(levelplot( m[,i] ~ plon+plat, bs, aspect="iso"))
      }

      # more range checks
      s = snowcrab_lbm( p=p, DS="baseline", ret="sd", varnames=varnames )
      # range checks
      s = log( exp(s[[1]]) * s[[2]] ) # s[[2]] is serving as weight/probabilities
      
      sq = quantile(s, probs=p$lbm_quantile_bounds[2], na.rm=TRUE ) 
      # s[which(s > sq)] = sq  # cap upper bound of sd

      # mm = which(( m - 1.96*s ) < 0 )
      # if (length(mm)>0) {
      #   m[mm] = NA
      #   s[mm] = NA 
      # }

      # limit range of extrapolation to within a given distance from survey stations .. annual basis
      set = snowcrab.db( DS="set.clean")
      bs = bathymetry.db(p=p, DS="baseline")
      bb = array_map( "xy->1", bs, gridparams=p$gridparams )

      if (0) {
        # annual mask 
        for (iy in 1:p$ny ) {
          S = set[ which(set$yr==p$yrs[iy]), c("plon", "plat") ] 
          S = S[ !duplicated(S),]
          nn = array_map( "xy->1", S, gridparams=p$gridparams )
          overlap = match(  nn, bb )
          overlap = overlap[ which( is.finite( overlap ))] 
          o = bs[overlap,] 
          # add corners as buffer
          ot = t(o) # reshape to make addition simpler using R's cycling rules
          o = rbind( t( ot + p$threshold.distance*c( 1, 1) ), 
                     t( ot + p$threshold.distance*c(-1,-1) ),
                     t( ot + p$threshold.distance*c( 1,-1) ),
                     t( ot + p$threshold.distance*c(-1, 1) )
          )
          o = o[ !duplicated(o),]
          boundary= non_convex_hull( o, alpha=p$threshold.distance*4, plot=FALSE )
          outside.polygon = which( point.in.polygon( bs[,1], bs[,2], boundary[,1], boundary[,2] ) == 0 )
          m[outside.polygon,iy] = NA
          s[outside.polygon,iy] = NA
          h[outside.polygon,iy] = NA
        }
      }

      # a static mask
      S = set[ , c("plon", "plat") ] 
      S = S[ !duplicated(S),]
      nn = array_map( "xy->1", S, gridparams=p$gridparams )
      overlap = match(  nn, bb )
      overlap = overlap[ which( is.finite( overlap ))] 
      o = bs[overlap,] 
      # add corners as buffer
      ot = t(o) # reshape to make addition simpler using R's cycling rules
      o = rbind( t( ot + p$threshold.distance*c( 1, 1) ), 
                 t( ot + p$threshold.distance*c(-1,-1) ),
                 t( ot + p$threshold.distance*c( 1,-1) ),
                 t( ot + p$threshold.distance*c(-1, 1) )
      )
      o = o[ !duplicated(o),]
      boundary= non_convex_hull( o, alpha=25, plot=FALSE )
      outside.polygon = which( point.in.polygon( bs[,1], bs[,2], boundary[,1], boundary[,2] ) == 0 )
      m[outside.polygon,] = NA
      s[outside.polygon,] = NA
      h[outside.polygon,] = NA

      B = list( m=m, s=s, h=h )
      save( B, file=fn, compress=TRUE )
      B = NULL

      projectdir = file.path(p$project.root, "maps", "fishable.biomass", p$spatial.domain )
      datarange = seq( qs[1]-0.1, qs[2]+0.1, length.out=150)
      cols = color.code( "seis", datarange )
      m[which(!is.finite(m))] = qs[1]
      for (iy in 1:p$ny) {
        y = p$yrs[iy]
        outfn = paste( "prediction.abundance.mean", y, sep=".")
        xyz = cbind( bs[, c("plon", "plat")], m[,iy] )
        map( xyz=xyz, cfa.regions=T, depthcontours=T, pts=NULL, annot=y,
          annot.cex=p$annot.cex, corners=p$planar.corners, fn=outfn, loc=projectdir, at=datarange,
          col.regions=cols, rez=c(p$pres,p$pres) )
      }

      datarange = seq( 0, max(s,na.rm=TRUE), length.out=150)
      cols = color.code( "seis", datarange )
      s[which(!is.finite(s))] = 0.1
      s[which(s < 0.1)] = 0.1
      
      for (iy in 1:p$ny) {
        y = p$yrs[iy]
        outfn = paste( "prediction.abundance.sd", y, sep=".")
        xyz = cbind( bs[, c("plon", "plat")], s[,iy] )
        map( xyz=xyz, cfa.regions=T, depthcontours=T, pts=NULL, annot=y,
          annot.cex=p$annot.cex, corners=p$planar.corners, fn=outfn, loc=projectdir, at=datarange,
          col.regions=cols, rez=c(p$pres,p$pres) )
      }

      return(fn)
    }

  
    # ------------------


    if (DS=="timeseries") {
      bm = interpolation.db(p=p, DS="biomass")
      bm$m = exp(bm$m) / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
      bm$s = exp(bm$s) / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
      
      bs = bathymetry.db( p=p, DS="baseline")
      
      bq = min( bm$m[ bm$m > 0 ], na.rm=T )

      K = NULL
      nreg = length(p$regions)
      for (r in 1:nreg ){
        aoi = filter.region.polygon(x=bs[ , c("plon", "plat")], region=p$regions[r], planar=T)
        aoi = intersect( aoi, which( bs$plon > 250 ) )
        out = matrix( NA, nrow=p$ny, ncol=3) 
        
        for (y in 1:p$ny) {
          iHabitat = which( bm$h[,y] > p$habitat.threshold.quantile ) # any area with biomass > lowest threshold
          iHabitatRegion = intersect( aoi, iHabitat )
          out[ y, 1] = sum( bm$m[iHabitatRegion,y] , na.rm=TRUE ) # abundance weighted by Pr
          out[ y, 2] = sqrt( sum( (bm$s[iHabitatRegion,y])^2 , na.rm=TRUE ) )
          out[ y, 3] = sum( bm$h[iHabitatRegion,y] ) * (p$pres*p$pres)
        }
        
        ok = as.data.frame( out )
        names( ok) = c("total", "total.sd", "sa.region")
        
        ok$log.total = log(ok$total)
        ok$total.sd.ln = log(ok$total.sd) # as above
        ok$region = p$regions[r]
        ok$yr = p$yrs

        ok$lbound = ok$total - ok$total.sd*1.96 # normal assumption 
        ok$ubound = ok$total + ok$total.sd*1.96 
        K = rbind(K, ok)
      }

      return( K )      
    }


    if ( DS %in% c( "interpolation.simulation" ) ) {
      message(" simulation-based results are not ready at present")
      message(" defaulting to simple estimates based upon assymptotic assumptions" )
      message(" only R0.mass is supported for now" )
      
      out = NULL  
      if ( p$vars.to.model == "R0.mass" ) out = interpolation.db( p=p, DS="timeseries" )
      
      return(out)
    }
 

  #############################

  # ---------


    if (DS =="habitat.temperatures") {

      bm = bio.snowcrab::interpolation.db( p=p, DS="biomass" )
      ps = bio.snowcrab::snowcrab_lbm(p=p, DS="output_data", voi=p$selection$name )

      temp = ps$t * bm$h

      K = NULL
      nreg = length(p$regions)
      for (r in 1:nreg ){
        aoi = filter.region.polygon(x=bs[ , c("plon", "plat")], region=p$regions[r], planar=T)
        aoi = intersect( aoi, which( bs$plon > 250 ) )
        out = matrix( NA, nrow=p$ny, ncol=2) 
        
        for (y in 1:p$ny) {
          iHabitat = which( bm$h[,y] > p$habitat.threshold.quantile ) # any area with bm > lowest threshold
          iHabitatRegion = intersect( aoi, iHabitat )
          out[ y, 1] = mean( temp[iHabitatRegion,y] , na.rm=TRUE ) # temperature weighted by Pr
          out[ y, 2] = sd( temp[iHabitatRegion,y] , na.rm=TRUE ) # temperature weighted by Pr
        }
        
        ok = as.data.frame( out )
        names( ok) = c("temperature", "temperature.sd")
        ok$region = p$regions[r]
        ok$yr = p$yrs
        ok$lbound = ok$temperature - ok$temperature.sd*1.96 # normal assumption 
        ok$ubound = ok$temperature + ok$temperature.sd*1.96 
        K = rbind(K, ok)
      }

      return( K )      

    }



    if ( DS %in% c( "interpolation.redo.old", "interpolation.old", "interpolation.simulation.old", "interpolation.simulation.redo.old",
                    "interpolation.simulation.complete.old", "interpolation.simulation.PS.old" ) ) {
    
  # from older methods .. keep for now till placed in a simulation method for lbm__habitat
   message( "deprecated .. for reference only until lbm__habitat is complete")

      basedir = file.path( project.datadirectory("bio.snowcrab"), "R", "gam" )

      loc.map = file.path( basedir, "maps" )
      loc.sol = file.path( basedir, "predictions" )
      loc.res = file.path( basedir, "results" )

      dir.create(path=loc.map, recursive=T, showWarnings=F)
      dir.create(path=loc.sol, recursive=T, showWarnings=F)
      dir.create(path=loc.res, recursive=T, showWarnings=F)


      if( DS=="interpolation.simulation.complete") {
        load( p$ofname )
        return( K )
      }

      if (DS=="interpolation.simulation.PS" ) {
        out = NULL
          v = p$v
          y = p$y
          fn.PS = file.path( loc.sol, paste( "PS.simulation.means", v, y, "rdata", sep="." ) )
          if ( ! (file.exists( fn.PS)) ) return(NULL)
          load(fn.PS)
        return (PS)
      }

      if (exists( "libs", p)) RLibrary( p$libs )
      if (is.null(ip)) ip = 1:p$nruns

      if (DS %in% c("interpolation.simulation")  ) {
        out = NULL
        for ( iip in ip ) {
          y = p$runs[iip,"y"]
          v = p$runs[iip,"v"]
          fn.K = file.path( loc.res, paste( "K", v, y, "rdata", sep="." ) )
          if ( ! (file.exists( fn.K)) ) next()
          load(fn.K)

				# this is a temporary fix for when data are not all properly refreshed
          oo = names(out)
					kk = names(K)

					if ( length(oo) > length(kk) ) {
						a = out
						b = K
					} else {
						a = K
						b = out
					}

					oo = setdiff( names(a), names(b) )
					if ( !is.null(out) && (length( oo ) > 0) ) {
						newvars = setdiff( names(a), names(b) )
						b[,newvars] = NA
						b = b[, names(a) ]
						out = rbind( a, b)

					} else {
						out = rbind( out, K )
					}

        }
        K = out
        save( K, file=p$ofname, compress=T )  # glue all the data results together into one file
        return ( K )
      }

      K = NULL

      for ( iip in ip ) {

        y = p$runs[iip,"y"]
        v = p$runs[iip,"v"]

        print ( p$runs[iip,] )

        qs = empirical.ranges( db="snowcrab", v, probs=c(p$habitat.threshold.quantile, 0.95), remove.zeros=T )

        PS = indicators.db ( DS="complete", year=y, p=p )
				PS$dyear = p$prediction.dyear  # must be same as above

        PST = temperature.db( p=p, DS="spacetime.prediction", yr=y  )
				if (is.null(PST)) next ()

        dyears = (c(1:(p$nw+1))-1)  / p$nw # intervals of decimal years... fractional year breaks
        pred.dyear.int = as.numeric( cut(p$prediction.dyear, breaks=dyears, include.lowest=T, ordered_result=TRUE ) ) # integerr representation of season
				PS$t = as.vector( PST[, pred.dyear.int ] )

        rm (PST); gc()

        PS$t[ which(PS$t < -2) ] = -2
			  PS$t[ which(PS$t > 30) ] = 30

        iitna = which( ! is.finite( PS$t ) )
        if (length(iitna)>0) PS$t[iitna] = PS$tmean[iitna]

        PS$z = log(PS$z)
        PS$dt.seasonal = PS$tmean - PS$t
        PS$dt.annual = PS$tmean - PS$tmean.climatology
        PS$sa = 1
        PS$dyear = p$prediction.dyear

        # PS$dZ = log( PS$dZ)
        # PS$ddZ = log( PS$ddZ)
        if (exists("tamp", PS)) PS$tamp = log(PS$tamp)
        if (exists("tamplitude.climatology", PS)) PS$tamplitude.climatology = log(PS$tamplitude.climatology)

				# posterior simulations
        Hmodel = NULL
        Hmodel = habitat.model.db( DS="habitat", v=v, yr=y )
        if (is.null( Hmodel)) next()

        Hsim = gam.simulation( M=Hmodel, X= PS, nsims=p$nsims ) #~8min
        rm( Hmodel); gc()

        print("finished H sim")

        oops = which( is.na(Hsim) )
        if (length(oops) > 0)  Hsim[oops ] = 0  # assume to be zero

        Amodel = NULL
			  Amodel = habitat.model.db( DS="abundance", v=v, yr=y )
        if (is.null(Amodel)) next()

        Asim = gam.simulation( M=Amodel, X= PS, nsims=p$nsims ) # ~5min
	      rm( Amodel); gc()
        print("finished A sim")

        oops = which( is.na(Asim) )
        if (length(oops) > 0)  Asim[oops ] = 0  # assume to be zero

        # Do not extrapolate: trim to XX% quantiles to be a little more conservative
        oopu =  which( Asim > qs[2] )
        if (length(oopu) > 0)  Asim[ oopu ] = qs[2]

        oopl =  which( Asim < qs[1]  )
        if (length(oopl) > 0)  Asim[ oopl ] = 0  # below detection limits

        Asim = Asim * Hsim  # Asim now becomes weighted by Pr of habitat

				Hsim.sa = colSums( Hsim ) # Pr weighted sum of habitat
				totalsurfacearea = mean ( Hsim.sa ) * (p$pres*p$pres)
        totalsurfacearea.sd = sd( Hsim.sa ) * (p$pres*p$pres)
        rm ( Hsim.sa ); gc()

        PS$habitat.mean = apply( Hsim, 1, mean, na.rm=T )
        PS$habitat.sd = apply( Hsim, 1, sd, na.rm=T )

        PS$abundance.mean = apply( Asim, 1, mean, na.rm=T )
        PS$abundance.sd =  apply( Asim, 1, sd, na.rm=T )

        # iAbundance = which ( PS$abundance.mean >= qs[1] )  #  consider abundance only if it is sufficiently precise (low associated variance)
        iHabitat = which( PS$habitat.mean > p$habitat.threshold.quantile  & (PS$habitat.mean - 2 * PS$habitat.sd) > 0 )
        iStations = filter.prediction.locations( DS="limit.to.near.survey.stations", PS=PS, y=y, p=p )  # consider locations only if close to a survey location (minimal extrapolation)
        # iHA = intersect(iHabitat, iAbundance)  # good locations for habitat and abundance prediction
				iE = union( iStations, intersect(iHabitat, iStations ))  #  good locations for habitat and abundance prediction, but filtered for proximity to survey stations
				rm( iStations); gc()

        # levelplot( habitat.mean ~ plon + plat, data=PS[,], aspect="iso" )
        # levelplot( habitat.mean ~ plon + plat, data=PS[iE,], aspect="iso" )

        fn.res = file.path( loc.sol, paste( "PS.simulation.means", v, y, "rdata", sep="." ) )
        print (fn.res )
        save( PS, file=fn.res, compress=T )

        K = NULL

        fn.K = file.path( loc.res, paste( "K", v, y, "rdata", sep="." ) )
        print(fn.K)

        PS$abundance.mean.log = log10( PS$abundance.mean )  # only used for plotting
        er = log10( qs  )

        PS$abundance.mean.log [ which( PS$abundance.mean.log < er[1] ) ] = er[1]
        # levelplot( abundance.mean.log ~ plon + plat, data=PS, aspect="iso" )

        if ("map.habitat" %in% p$ItemsToMap ) {
          datarange = seq( 0, 1, length.out=150)
          cols = color.code( "seis", datarange )
          outfn = paste( "prediction.habitat.mean", v, y, sep=".")
          map( xyz=PS[,c("plon", "plat", "habitat.mean")], cfa.regions=T, depthcontours=T, pts=NULL, annot=paste( v, y ),
            annot.cex=p$annot.cex, corners=p$planar.corners, fn=outfn, loc=loc.map, at=datarange,
            col.regions=cols, rez=c(p$pres,p$pres) )
        }

        if ("map.abundance" %in% p$ItemsToMap ) {
          datarange = seq( er[1], er[2], length.out=150)
          cols = color.code( "seis", datarange )
          outfn = paste( "prediction.abundance.mean", v, y, sep=".")
          map( xyz=PS[ , c("plon", "plat", "abundance.mean.log")], cfa.regions=T, depthcontours=T, pts=NULL, annot= paste( v, y ),
            annot.cex=p$annot.cex, corners=p$planar.corners, fn=outfn, loc=loc.map, at=datarange,
            col.regions=cols, rez=c(p$pres,p$pres) )
        }

        if ("map.abundance.estimation" %in% p$ItemsToMap ) {
          tomap = NULL
          tomap = PS[ , c("plon", "plat", "abundance.mean.log")]
          tomap$abundance.mean.log[ setdiff( 1:nrow(PS), iE ) ] = log10( er[1] ) - 1  # just a value smaller than the lower bound
          datarange = seq( er[1], er[2], length.out=150)
          cols = color.code( "seis", datarange )
          outfn = paste( "prediction.abundance.mean.estimationarea", v, y, sep=".")
          map( xyz=tomap, cfa.regions=T, depthcontours=T, pts=NULL,
            annot=paste( v, y ), annot.cex=p$annot.cex, corners=p$planar.corners, fn=outfn, loc=loc.map, at=datarange,
            col.regions=cols, rez=c(p$pres,p$pres) )
          rm (tomap) ; gc()
        }

        PS$abundance.mean.log = NULL
				gc()

        print( "finished maps")

        for (r in p$regions ){

          iRegion.PS = filter.region.polygon(x=PS[ , c("plon", "plat")], region=r, planar=T)
          iEastof250 = which( PS$plon > 250 )
          iRegion.PS = intersect( iRegion.PS, iEastof250 )

					iHabitatSubarea = intersect( iRegion.PS, iHabitat ) # plotting surface and real habitat area
          iEstimationArea = intersect( iRegion.PS, iE )

          if ( length( iEstimationArea ) > 10 ) {

            hhsum = colSums(  Hsim[ iHabitatSubarea , ] ) # area weighted by Pr
            sa.region = mean( hhsum ) * (p$pres*p$pres)
            sa.sd = sd( hhsum ) * (p$pres*p$pres)

            hhsum = colSums( Hsim[ iEstimationArea , ] )
            sa.estim = mean( hhsum ) * (p$pres*p$pres)
            sa.estim.sd = sd( hhsum ) * (p$pres*p$pres)

            aa.sums = apply( Asim[ iEstimationArea,  ] , 2, sum ) # abundance weighted by Pr
            V.mean = mean( aa.sums )
            V.sd = sd( aa.sums )
            V.sd.ln = sd( log( aa.sums ), na.rm=T )
            ci = quantile( aa.sums, probs=c(0.025, 0.5, 0.975), na.rm=T,names=F )

          } else {
            V.mean = NA
            V.sd = NA
            V.sd.ln = NA
            sa.estim = NA
            sa.estim.sd = NA
            sa.region = NA
            sa.sd = NA
            ci = rep( NA, 3 )
          }

          bb.sums = apply( Asim[ iHabitatSubarea,  ] , 2, sum ) # abundance weighted by Pr
          W.mean = mean( bb.sums )
          W.sd = sd( bb.sums )
          W.ci = quantile( bb.sums, probs=c(0.025, 0.5, 0.975), na.rm=T,names=F )

          L = data.frame( yr=y, vars=v, region=r,
              total=V.mean, total.sd=V.sd, total.sd.ln=V.sd.ln, median=ci[2], lbound=ci[1], ubound=ci[3],
              ss=W.mean, ss.sd=W.sd, ss.median=W.ci[2], ss.lbound=W.ci[1], ss.ubound=W.ci[3],
              sa.total=totalsurfacearea, sa.total.sd=totalsurfacearea.sd,
              sa.region=sa.region, sa.region.sd=sa.sd,
              sa.estimation=sa.estim, sa.estimation.sd=sa.estim.sd,

              datestamp=as.character( Sys.time() ) )
          K = rbind( K, L )
        } # regions

        save(K, file=fn.K, compress=T)
        print (K)
        print( "... Completed successfully" )

      } # runs (ip)
    } # end if interpolations...

    return ("Completed" )
  }


