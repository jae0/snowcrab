
snowcrab.habitat.db = function( ip=NULL, DS=NULL, p=NULL, voi=NULL, y=NULL, selection=NULL ) {

  # over-ride default dependent variable name if it exists
  if (is.null(voi)) if (!is.null(selection)) if (exists("name", selection)) voi=selection$name
  if (is.null(voi)) if (exists("variables",p)) if(exists("Y", p$variables)) voi=p$variables$Y

  if (DS %in% c("baseline", "baseline.redo") ) {
    # based upon bio.snowcrab::habitat.model.db()

    outdir = file.path( project.datadirectory("bio.snowcrab"), "lbm"  )
    dir.create(path=outdir, recursive=T, showWarnings=F)
    fn = file.path( outdir, paste("baseline", voi, "rdata", sep=".") )

    set = NULL
    if ( DS == "baseline" ) {
      if (file.exists(fn)) load(fn)
      return( set)
    }

    set = bio.indicators::survey.db( p=p, DS="set.filter", selection=selection ) # mature male > 74 mm 

    set = presence.absence( X=set, vname="zm", px=p$habitat.threshold.quantile )  # determine presence absence and weighting

    # add commerical fishery data
    lgbk = logbook.db( DS="fisheries.complete", p=p )
    lgbk = lgbk[ which( is.finite( lgbk$landings)), ]
    lgbk$totmass = NA # dummy to bring in mass as well 
    lgbk$data.source = "logbooks"
    lgbk$lon = lgbk$lat = NULL
        
    lgbk = presence.absence( X=lgbk, vname="landings", px=p$habitat.threshold.quantile )  # determine presence absence and weighting

    # baddata = which( lgbk$z < log(50) | lgbk$z > log(600) )
    # if ( length(baddata) > 0 ) lgbk = lgbk[ -baddata,]

    # lgbk$julian = lubridate::yday( lgbk$date.landed )

    nms = intersect( names(set) , names( lgbk) )
    set = rbind( set[, nms], lgbk[,nms] )

    names(set)[ which( names(set) =="totmass")] = selection$name 


    # Z = bathymetry.db( DS="baseline", p=p )
    # Z$plon = floor(Z$plon / 10)*10
    # Z$plat = floor(Z$plat / 10)*10
    # ii = which(duplicated(Z))
    # if (length(ii)>0) Z = Z[-ii,] # thinned list of locations

    # dd = rdist( set[,c("plon", "plat")] , Z )
    # ee = apply( dd, 1, min, na.rm=T )
    # ff = which( ee < p$threshold.distance ) # all within XX km of a good data point
    # set = set[ ff, ]

    # # bring in time invariant features:: depth
    # print ("Bring in depth")
    # set = habitat.lookup( set,  p=p, DS="depth" )
    # set$z = log( set$z )

    # # bring in time varing features:: temperature
    # print ("Bring in temperature")
    # set = habitat.lookup( set, p=p, DS="temperature" )

    # # bring in all other habitat variables, use "z" as a proxy of data availability
    # # and then rename a few vars to prevent name conflicts
    # set = habitat.lookup( set,  p=p, DS="all.data" )

    # return planar coords to correct resolution
    # set = lonlat2planar( set, proj.type=p$internal.projection )

    # # complete area designations
    # set = fishing.area.designations(set, type="lonlat")



    save ( set, file=fn, compress=TRUE )

    return (fn)

  }



  if (DS=="lbm_inputs") {
    # mostly based on indicators.db( DS="lbm_inputs") 

    INP = snowcrab.habitat.db(p=p, DS="baseline", voi=selection$name )
    INP$tiyr = lubridate::decimal_date( INP$timestamp ) 

    locsmap = match( 
      lbm::array_map( "xy->1", INP[,c("plon","plat")], gridparams=p$gridparams ), 
      lbm::array_map( "xy->1", bathymetry.db(p=p, DS="baseline"), gridparams=p$gridparams ) )

    # spatial vars and climatologies 
    newvars = c("dZ", "ddZ", "log.substrate.grainsize", "tmean", "tsd" )
    sn = indicators.lookup( p=p, DS="spatial", locsmap=locsmap, varnames=newvars )
    names( sn ) = c("dZ", "ddZ", "log.substrate.grainsize", "tmean.climatology", "tsd.climatology" )
    INP = cbind( INP,  sn )

    # for space-time(year-averages) 
    newvars = c( "tmean", "tsd", "amplitude" )
    sn = indicators.lookup( p=p, DS="spatial.annual", locsmap=locsmap, timestamp=INP[,"timestamp"], varnames=newvars )
    colnames( sn  ) = newvars
    INP = cbind( INP,  sn )
    INP$tamplitude = INP$amplitude

    # additional indicators.db variables
    for (iv in names(p$indicators.variables)) {
      p0 = bio.indicators::indicators.parameters( p=p, DS="default", current.year=p$current.year )
      p0 = bio.indicators::indicators.parameters( p=p0, DS=iv  )
      p0 = bio.spacetime::spatial_parameters( p=p0, type=p$spatial.domain ) # return to correct domain
      vn = p0$indicators.variables[iv]
      sn = indicators.lookup( p=p0, DS="spatial.annual", locsmap=locsmap, timestamp=INP[,"timestamp"], 
        varnames=vn, DB=indicators.db( p=p0, DS="baseline", varnames=vn ) )
      colnames( sn  ) = p$indicators.variables[iv]
      INP = cbind( INP,  sn )
    }

    INP = INP[, which(names(INP) %in% c(p$varnames, p$variables$Y, p$variables$TIME, "dyear", "yr" ) ) ]  # a data frame
    oo = setdiff(p$varnames, names(INP))
    if (length(oo) > 0 ) {
      print(oo )
      warning("Some variables are missing in the input data")
    }
    INP = na.omit(INP)

    # cap quantiles of dependent vars
      dr = list()
      for (voi in p$varnames) {
        dr[[voi]] = quantile( INP[,voi], probs=p$lbm_quantile_bounds, na.rm=TRUE ) # use 95%CI
        il = which( INP[,voi] < dr[[voi]][1] )
        if ( length(il) > 0 ) INP[il,voi] = dr[[voi]][1]
        iu = which( INP[,voi] > dr[[voi]][2] )
        if ( length(iu) > 0 ) INP[iu,voi] = dr[[voi]][2]
      }

      PS = indicators.db( p=p, DS="prediction.surface" ) # a list object with static and annually varying variables  
      names(PS)[ names(PS)=="amplitude"] ="tamplitude" 

      ps_varnames = setdiff( p$varnames, p$variables$LOCS )

      PS = PS[ which(names(PS) %in% ps_varnames ) ] # time vars, if they are part of the model will be created within lbm

      oo = setdiff(p$varnames, ps_varnames )
      if (length(oo) > 0 ) {
        print(oo )
        warning("Some variables are missing in the prediction surface, PS")
      }

      OUT = list( LOCS=bathymetry.db(p=p, DS="baseline"), COV=PS )    

      return (list(input=INP, output=OUT))

  }

  
  # ----------------

  
  if (DS %in% c("prediction.surface", "prediction.surface.redo") ) {

    outdir = file.path( project.datadirectory("bio.snowcrab"), "PS", p$spatial.domain )
    dir.create(outdir, recursive=T, showWarnings=F)

    dyear_index = 1
    if (exists("dyears", p) & exists("prediction.dyear", p))  dyear_index = which.min( abs( p$prediction.dyear - p$dyears))

    outfile =  file.path( outdir, paste("PS", dyear_index, "rdata", sep=".") )

    if ( DS=="prediction.surface" ) {
      PS = NULL
      if (file.exists(outfile)) load( outfile )
      return (PS)
    }

    # this is the same as the p`rediction surface as used for indicators.db but the domain is smaller
    # .. more will be added to it belwo
    PS = indicators.db(p=p, DS="spatial")
    names(PS)[which(names(PS)=="tmean")] = "tmean.climatology"
    names(PS)[which(names(PS)=="tsd")] = "tsd.climatology"
    names(PS)[which(names(PS)=="tmin")] = "tmin.climatology"
    names(PS)[which(names(PS)=="tmax")] = "tmax.climatology"
    names(PS)[which(names(PS)=="amplitude")] = "tamplitude.climatology"

    nPS = nrow( PS )
    PS = as.list(PS)

    p0 = bio.temperature::temperature.parameters(p=p, current.year=p$current.year )
    p0 = bio.temperature::temperature.parameters( DS="lbm", p=p0 )
    p0 = bio.spacetime::spatial_parameters( p=p0, type=p$spatial.domain ) # return to correct domain

    yr_index = match( p$yrs, p0$tyears )
    u = indicators.db(p=p, DS="spatial.annual")
    for ( vn in names(u) ){
      u[[vn]] = u[[vn]][,yr_index]
    }
    PS = c( PS, u)

    # now we add the other covariate fields for modelling and prediction
    # additional indicators.db variables
    for (iv in names(p$indicators.variables)) {
      p0 = bio.indicators::indicators.parameters( p=p, DS="default", current.year=p$current.year )
      p0 = bio.indicators::indicators.parameters( p=p0, DS=iv  )
      p0 = bio.spacetime::spatial_parameters( p=p0, type=p$spatial.domain ) # return to correct domain
      DB=indicators.db( p=p0, DS="baseline", varnames=p0$indicators.variables[iv] )
      yr_index = match( p$yrs, p0$yrs )
      PS = c( PS, DB[,yr_index] )
    }

    save (PS, file=outfile, compress=T )
    return( outfile )

  }




  outdir = file.path( project.datadirectory("bio.snowcrab"), "R", "gam", "habitat" )
  dir.create(path=outdir, recursive=T, showWarnings=F)

  if (DS=="PS" ) {
    out = NULL
      fn.PS = file.path( outdir, paste( "PS", voi, y, "rdata", sep="." ) )
      if ( ! (file.exists( fn.PS)) ) return(NULL)
      load(fn.PS)
    return (PS)
  }

  if (DS=="K" ) {
    out = NULL
      fn.K = file.path( outdir, paste( "K", voi, y, "rdata", sep="." ) )
      if ( ! (file.exists( fn.K)) ) return(NULL)
      load(fn.K)
    return (K)
  }

  if (!is.null(p$libs)) for( i in p$libs ) require(i)
  if (is.null(ip)) ip = 1:p$nruns

  for ( iip in ip ) {

      y = p$runs[iip,"y"]
      v = p$runs[iip,"v"]

      v0 = v = p$runs[iip,"v"]
      if ( v0 =="R0.mass.environmentals.only" ) v="R0.mass"


      PS = indicators.db ( DS="complete", year=y, p=p )
			PS$dyear = p$prediction.dyear  # must be same as above
			PS$t = NA

      PST = temperature.db( p=p, DS="spacetime.prediction", yr=y  )
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
      PS$dt.annual = PS$tmean - PS$tmean.climatology
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
          annot.cex=p$annot.cex, corners=p$planar.corners, fn=outfn, loc=outdir, at=datarange,
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



