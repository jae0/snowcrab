
snowcrab_lbm = function( ip=NULL, DS=NULL, p=NULL, voi=NULL, y=NULL, selection=NULL ) {

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

  # -----------------------------------

  if (DS=="lbm_inputs") {
    # mostly based on indicators.db( DS="lbm_inputs") 

    INP = snowcrab_lbm(p=p, DS="baseline", voi=selection$name )
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
   # mostly based on indicators.db( DS="prediction.surface") 

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


  # --------------------------


    if ( DS %in% c("predictions", "predictions.redo" ) ) {
      # NOTE: the primary interpolated data were already created by lbm. 
      # This routine points to this data and also creates 
      # subsets of the data where required, determined by "spatial.domain.subareas" 
      
      projectdir = file.path(p$project.root, "modelled", voi, p$spatial.domain )
      
      if (DS %in% c("predictions")) {
        P = V = NULL
        fn = file.path( projectdir, paste("lbm.prediction", ret,  year, "rdata", sep=".") )
        if (is.null(ret)) ret="mean"
        if (file.exists(fn) ) load(fn) 
        if (ret=="mean") return (P)
        if (ret=="sd") return( V)
      }

      if (exists( "libs", p)) RLibrary( p$libs )
      if ( is.null(ip) ) ip = 1:p$nruns

      # downscale and warp from p(0) -> p1

      for ( r in ip ) {
        # print (r)
        year = p$runs[r, "yrs"]
        # default domain
        PP0 = lbm_db( p=p, DS="lbm.prediction", yr=year, ret="mean")
        VV0 = lbm_db( p=p, DS="lbm.prediction", yr=year, ret="sd")
        p0 = spatial_parameters( p=p, type=p$spatial.domain ) # from
        L0 = bathymetry.db( p=p0, DS="baseline" )
        L0i = lbm::array_map( "xy->2", L0[, c("plon", "plat")], gridparams=p0$gridparams )
        sreg = setdiff( p$spatial.domain.subareas, p$spatial.domain ) 

        for ( gr in sreg ) {
          p1 = spatial_parameters( p=p, type=gr ) # 'warping' from p -> p1
          L1 = bathymetry.db( p=p1, DS="baseline" )
          L1i = lbm::array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
          L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
          L1$plon_1 = L1$plon # store original coords
          L1$plat_1 = L1$plat
          L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
          p1$wght = fields::setup.image.smooth( nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres, theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
          P = spatial_warp( PP0[], L0, L1, p0, p1, "fast", L0i, L1i )
          V = spatial_warp( VV0[], L0, L1, p0, p1, "fast", L0i, L1i )
          projectdir_p1 = file.path(p$project.root, "modelled", voi, p1$spatial.domain ) 
          dir.create( projectdir_p1, recursive=T, showWarnings=F )
          fn1_sg = file.path( projectdir_p1, paste("lbm.prediction.mean",  year, "rdata", sep=".") )
          fn2_sg = file.path( projectdir_p1, paste("lbm.prediction.sd",  year, "rdata", sep=".") )
          save( P, file=fn1_sg, compress=T )
          save( V, file=fn2_sg, compress=T )
          print (fn1_sg)
        }
      } 
      return ("Completed")

      if (0) {
        levelplot( P ~ plon_1 + plat_1, L1, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      }
    
    }



    #  -------------------------------

    if (DS %in% c(  "lbm.stats", "lbm.stats.redo" )){

      
      if (DS %in% c("lbm.stats")) {
        stats = NULL
        projectdir = file.path(p$project.root, "modelled", voi, p$spatial.domain )
        fn = file.path( projectdir, paste( "lbm.statistics", "rdata", sep=".") )
        if (file.exists(fn) ) load(fn) 
        return( stats )
      }

      # downscale and warp from p(0) -> p1
      # default domain
      S0 = lbm_db( p=p, DS="stats.to.prediction.grid" )
      Snames = colnames(S0)
      p0 = spatial_parameters( p=p, type=p$spatial.domain ) # from
      L0 = bathymetry.db( p=p0, DS="baseline" )
      L0i = lbm::array_map( "xy->2", L0[, c("plon", "plat")], gridparams=p0$gridparams )
      sreg = setdiff( p$spatial.domain.subareas, p$spatial.domain ) 

      for ( gr in sreg ) {
        p1 = spatial_parameters( p=p, type=gr ) # 'warping' from p -> p1
        L1 = bathymetry.db( p=p1, DS="baseline" )
        L1i = lbm::array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
        L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
        L1$plon_1 = L1$plon # store original coords
        L1$plat_1 = L1$plat
        L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
        p1$wght = fields::setup.image.smooth( nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres, theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
        stats = matrix( NA, ncol=ncol(S0), nrow=nrow(L1) )
        for ( i in 1:ncol(S0) ) {
          stats[,i] = spatial_warp( S0[,i], L0, L1, p0, p1, "fast", L0i, L1i )
        }
        colnames(stats) = Snames
        projectdir_p1 = file.path(p$project.root, "modelled", voi, p1$spatial.domain ) 
        dir.create( projectdir_p1, recursive=T, showWarnings=F )
        fn1_sg = file.path( projectdir_p1, paste("lbm.statistics", "rdata", sep=".") )
        save( stats, file=fn1_sg, compress=T )
        print (fn1_sg)
      }
      return ("Completed")
    
      if (0) {
        levelplot( stats[,1] ~ plon_1 + plat_1, L1, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      }
   

    }


    #  -------------------------------


    if (DS %in% c("complete", "complete.redo") ) {
      # assemble data for a given project 
   
      if (DS=="complete") {
        IC = NULL
        projectdir = file.path(p$project.root, "modelled", voi, p$spatial.domain )
        dir.create(projectdir, recursive=T, showWarnings=F)
        outfile =  file.path( projectdir, paste( "snowcrab", "complete", p$spatial.domain, "rdata", sep= ".") )
        if ( file.exists( outfile ) ) load( outfile )
        Inames = names(IC)
        if (is.null(varnames)) varnames=Inames
        varnames = intersect( Inames, varnames )
        if (length(varnames) == 0) varnames=Inames  # no match .. send all
        IC = IC[ , varnames]
        return(IC)
      }

      if (exists( "libs", p)) RLibrary( p$libs )

      grids = unique( c(p$spatial.domain.subareas , p$spatial.domain ) ) # operate upon every domain
   
      for (gr in grids ) {
        print(gr)

        p1 = spatial_parameters( p=p, type=gr ) #target projection    
        L1 = bathymetry.db(p=p1, DS="baseline")

        BS = snowcrab_lbm( p=p1, DS="lbm.stats" )
        colnames(BS) = paste(voi, colnames(BS), sep=".")
        IC = cbind( L1, BS )

        # climatology
        nL1 = nrow(L1)
        PS = PSsd = matrix( NA, nrow=nL1, ncol=p$ny )
        if (is.null(voi)) p1$variables$Y = voi # need to send this to get the correct results
        for (iy in 1:p$ny) {
          yr = p$yrs[iy]
          PS[,iy] = lbm_db( p=p1, DS="lbm.prediction", yr=yr, ret="mean")
          PSsd[,iy] = lbm_db( p=p1, DS="lbm.prediction", yr=yr, ret="sd")
        }

        qPS = quantile( PS, probs=p$lbm_quantile_bounds, na.rm=TRUE )
        u = which( PS < qPS[1])
        if (length(u)>0) PS[u] = qPS[1]
        v = which( PS > qPS[2])
        if (length(v)>0) PS[v] = qPS[2]
        
        qPSsd = quantile( PSsd, probs=p$lbm_quantile_bounds, na.rm=TRUE )
        u = which( PSsd < qPSsd[1])
        if (length(u)>0) PSsd[u] = qPSsd[1]
        v = which( PSsd > qPSsd[2])
        if (length(v)>0) PSsd[v] = qPSsd[2]
      
        CL = cbind( apply( PS, 1, mean, na.rm=TRUE ), apply( PSsd, 1, mean, na.rm=TRUE ) )
        colnames(CL) = paste( voi, c("mean", "sd"), "climatology", sep=".")
        IC = cbind( IC, CL )
        PS = PSsd = NULL

        projectdir = file.path(p$project.root, "modelled", voi, p1$spatial.domain )
        dir.create( projectdir, recursive=T, showWarnings=F )
        outfile =  file.path( projectdir, paste( "snowcrab", "complete", p1$spatial.domain, "rdata", sep= ".") )
        save( IC, file=outfile, compress=T )
        print( outfile )

      }
      
      return( "Complete" )
    }

    # -------------------

    if (DS %in% c("baseline", "baseline.redo") ) {


      if ( DS=="baseline" ) {
        BL = list()
        for (voi in varnames ) {
          projectdir = file.path(p$project.root, "modelled", voi, p$spatial.domain )
          outfile =  file.path( projectdir, paste( "snowcrab", "baseline", p$spatial.domain, "rdata", sep= ".") )
          TS = NULL
          load( outfile)
          BL[[voi]] = TS
        }
        return (BL)
      }

      if (exists( "libs", p)) RLibrary( p$libs )
      if (is.null(ip)) ip = 1:p$nruns

      p$variables$Y = voi # need to send this to get the correct results
      grids = unique( c(p$spatial.domain.subareas , p$spatial.domain ) ) # operate upon every domain
   
      for (gr in grids ) {
        print(gr)
        p1 = spatial_parameters( p=p, type=gr ) #target projection    
        L1 = bathymetry.db(p=p1, DS="baseline")
        nL1 = nrow(L1)
        TS = matrix( NA, nrow=nL1, ncol=p$ny )
        for (i in 1:p$ny ) {
          yr = p$yrs[i]
          TS[,i] = lbm_db( p=p1, DS="lbm.prediction", yr=yr, ret="mean")
         }

        qTS = quantile( TS, probs=p$lbm_quantile_bounds, na.rm=TRUE )
        u = which( TS < qTS[1])
        if (length(u)>0) TS[u] = qTS[1]
        v = which( TS > qTS[2])
        if (length(v)>0) TS[v] = qTS[2]
  

        projectdir = file.path(p$project.root, "modelled", voi, p1$spatial.domain )
        dir.create( projectdir, recursive=T, showWarnings=F )
        outfile =  file.path( projectdir, paste( "snowcrab", "baseline", p1$spatial.domain, "rdata", sep= ".") )
        save( TS, file=outfile, compress=T )
        print( outfile )
      }

      return( "Complete" )
 
    }


  #############################


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



