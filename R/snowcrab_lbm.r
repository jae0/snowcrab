
snowcrab_lbm = function( ip=NULL, DS=NULL, p=NULL, voi=NULL, year=NULL, ret=NULL, varnames=NULL ) {

  # over-ride default dependent variable name if it exists
  if (is.null(voi)) if (exists("selection",p)) if (exists("name", p$selection)) voi=p$selection$name
  if (is.null(voi)) if (exists("variables",p)) if (exists("Y", p$variables))    voi=p$variables$Y
  if (exists( "libs", p)) RLibrary( p$libs )
   

  if (DS=="lbm_inputs") {
    # mostly based on indicators.db( DS="lbm_inputs") 
    INP = snowcrab_lbm(p=p, DS="input_data" )  # , voi=p$selection$name
    PS  = snowcrab_lbm(p=p, DS="output_data"  ) # , voi=p$selection$name
    LOCS = bathymetry.db(p=p, DS="baseline")
    return (list(input=INP, output=list( LOCS=LOCS, COV=PS )) )
  }


  # --------------------------

  if (DS %in% c("input_data") ) {
    set = bio.indicators::survey.db( p=p, DS="set.filter" ) # mature male > 95 mm 
 
    if ( p$selection$type=="abundance") {
      # snowcrab survey data only
      set = set[ which(set$data.source == "snowcrab"), ]
      
      # robustify input data: .. upper bound trim
      qq = quantile( set$totmass, probs=0.975, na.rm=TRUE )
      set$totmass[ set$totmass > qq] = qq

      # keep "zero's" to inform spatial processes but only as "lowestpossible" value
      jj = which( set$totmass > 0 )
      # set = set[jj,]
      lowestpossible = min( set$totmass[jj] , na.rm=TRUE)  
      ii = which( set$totmass < lowestpossible )
      set$totmass[ii] = lowestpossible / 4
      set$totmass = log( set$totmass)
      names(set)[ which( names(set) =="totmass")] = p$selection$name 
      set$Y = NULL
      set$wt = set$sa 
    }

    if ( p$selection$type=="presence_absence") {
      # must run here as we need the wgt from this for both PA and abundance
      set = presence.absence( X=set, vname="zm", px=p$habitat.threshold.quantile )  # determine presence absence and weighting

      if ( grepl( "snowcrab.large.males", p$selection$name ) ) {
        # add commerical fishery data -- 
        # depth data is problematic ... drop for now
        lgbk = logbook.db( DS="fisheries.complete", p=p )
        lgbk = lgbk[ which( is.finite( lgbk$landings)), ]
        lgbk$totmass = NA # dummy to bring in mass as well 
        lgbk$data.source = "logbooks"
        lgbk = presence.absence( X=lgbk, vname="landings", px=p$habitat.threshold.quantile )  # determine presence absence and weighting
        lgbk$z = exp( lgbk$z )
        nms = intersect( names(set) , names( lgbk) )
        set = rbind( set[, nms], lgbk[,nms] )
      }

      set = set[ which(set$data.source %in% c("snowcrab", "groundfish", "logbooks") ), ]
      names(set)[ which( names(set) =="Y")] = p$selection$name 
      set$totmass = NULL
      set = set[ which(is.finite(set$plon + set$plat)),]
    }

    set = set[ which(is.finite(set[, p$selection$name])),]

    set = set[ which(set$yr %in% p$yrs ), ]

    coast = bio.polygons::coastline.db( p=p, DS="mapdata.coastPolygon" )
    coast = spTransform( coast, CRS("+proj=longlat +datum=WGS84") )
    setcoord = SpatialPoints( as.matrix( set[, c("lon", "lat")]),  proj4string=CRS("+proj=longlat +datum=WGS84") )
    inside = sp::over( setcoord, coast )
    onland = which (is.finite(inside))
    if (length(onland)>0) set = set[-onland, ]
    
    set$lon = set$lat = NULL

    set$tiyr = lubridate::decimal_date( set$timestamp ) 

    bathy = bathymetry.db(p=p, DS="baseline")

    psse = p
    psse$spatial.domain = "SSE"
    psse = spatial_parameters( p=psse )  # data are from this domain .. so far
    bathysse = bathymetry.db(p=psse, DS="baseline")

    # redo as set has changed
    locsmap = match( 
      lbm::array_map( "xy->1", set[,c("plon","plat")], gridparams=p$gridparams ), 
      lbm::array_map( "xy->1", bathy[,c("plon","plat")], gridparams=p$gridparams ) )

    # spatial vars and climatologies 
    newvars = c("dZ", "ddZ", "log.substrate.grainsize", "tmean.climatology", "tsd.climatology", "b.range", "t.range" )
    sn = indicators.lookup( p=p, DS="spatial", locsmap=locsmap, varnames=newvars )
    set = cbind( set,  sn )

    oo = which( !is.finite(locsmap) )
    if (length(oo) > 0) {
      # catch stragglers from a larger domain 
      locsmapsse = match( 
        lbm::array_map( "xy->1", set[oo, c("plon","plat")], gridparams=psse$gridparams ), 
        lbm::array_map( "xy->1", bathysse[,c("plon","plat")], gridparams=psse$gridparams ) )
      sn = indicators.lookup( p=psse, DS="spatial", locsmap=locsmapsse, varnames=newvars )
      for (nv in newvars) set[oo,nv] = sn[,nv]
    }

    # for space-time(year-averages) 
    newvars = c( "tmean", "tsd", "amplitude" )
    sn = indicators.lookup( p=p, DS="spatial.annual", locsmap=locsmap, timestamp=set[,"timestamp"], varnames=newvars )
    colnames( sn  ) = newvars
    set = cbind( set,  sn )
    
    nn = which( !is.finite(set$tmean) )
    if (length(nn) > 0) {
      # catch stragglers from a larger domain 
      locsmapsse = match( 
        lbm::array_map( "xy->1", set[nn, c("plon","plat")], gridparams=psse$gridparams ), 
        lbm::array_map( "xy->1", bathysse[,c("plon","plat")], gridparams=psse$gridparams ) )
      sn = indicators.lookup( p=psse, DS="spatial.annual", locsmap=locsmapsse, timestamp=set[nn,"timestamp"], varnames=newvars )
      for (nv in newvars) set[nn,nv] = sn[,nv]
    }

    names(set)[ names(set)=="amplitude"] ="tamplitude"

    # additional indicators.db variables
    for (iv in names(p$indicators.variables)) {
      p0 = bio.indicators::indicators.parameters( p=p, DS="default", year.assessment=p$year.assessment )
      p0 = bio.indicators::indicators.parameters( p=p0, DS=iv  )
      p0 = lbm::spatial_parameters( p=p0, type=p$spatial.domain ) # return to correct domain
      vn = p0$indicators.variables[[iv]]
      sn = indicators.lookup( p=p0, DS="spatial.annual", locsmap=locsmap, timestamp=set[,"timestamp"], 
        varnames=vn, DB=indicators.db( p=p0, DS="baseline", varnames=vn ) )
      sn = as.data.frame(sn)
      names( sn  ) = p$indicators.variables[[iv]]
      set = cbind( set,  sn )

      mm = which( !is.finite(set[,vn[1]]) )
      if (length(mm) > 0) {
        # catch stragglers from a larger domain
        p0 = lbm::spatial_parameters( p=p0, type="SSE" ) 
        locsmapsse = match( 
          lbm::array_map( "xy->1", set[mm, c("plon","plat")], gridparams=p0$gridparams ), 
          lbm::array_map( "xy->1", bathysse[,c("plon","plat")], gridparams=p0$gridparams ) )
        sn = indicators.lookup( p=p0, DS="spatial.annual", locsmap=locsmapsse, timestamp=set[mm,"timestamp"], 
          varnames=vn, DB=indicators.db( p=p0, DS="baseline", varnames=vn ) )
        for (nv in vn) set[mm,nv] = sn[,nv]
      }
    }

    set = set[, which(names(set) %in% c(p$varnames, p$variables$Y, p$variables$TIME, "dyear", "yr",  "wt") ) ]  # a data frame
    oo = setdiff(p$varnames, names(set))
    if (length(oo) > 0 ) {
      print(oo )
      warning("Some variables are missing in the input data")
    }
    set = na.omit(set)

    # cap quantiles of dependent vars
    dr = list()
    ps_varnames = setdiff( p$varnames, p$variables$LOCS )
    
    for (pvn in ps_varnames) {
      dr[[pvn]] = quantile( set[,pvn], probs=p$lbm_quantile_bounds, na.rm=TRUE ) # use 95%CI
      il = which( set[,pvn] < dr[[pvn]][1] )
      if ( length(il) > 0 ) set[il,pvn] = dr[[pvn]][1]
      iu = which( set[,pvn] > dr[[pvn]][2] )
      if ( length(iu) > 0 ) set[iu,pvn] = dr[[pvn]][2]
    }

    return (set)

  }


  # ------------------------


  if (DS %in% c("output_data") ) {
    PS = indicators.db( p=p, DS="prediction.surface" ) # a list object with static and annually varying variables  
    names(PS)[ names(PS)=="amplitude"] ="tamplitude" 

    # make years coherent for temperatures
    p0 = bio.indicators::indicators.parameters(p=p, year.assessment=p$year.assessment )
    yr_index = match( p$yrs, p0$yrs )
    yg = which(is.finite(yr_index))
    ym = which(is.na(yr_index))
    if (length(ym) > 0) {
      for ( vn in c("t", p0$bstats) ) {
        PS[[vn]][yg] = PS[[vn]][,yr_index[yg]]
        PS[[vn]][ym] = rowMeans( PS[[vn]][], na.rm=TRUE )
      }
    } else {
      for ( vn in c("t", p0$bstats) ) {
        PS[[vn]][] = PS[[vn]][,yr_index]
      }
    }

    # indicators.db variables 
    for (iv in names(p$indicators.variables)) {
      p0 = bio.indicators::indicators.parameters( p=p, DS="default", year.assessment=p$year.assessment )
      p0 = bio.indicators::indicators.parameters( p=p0, DS=iv  )
      p0 = lbm::spatial_parameters( p=p0, type=p$spatial.domain ) # return to correct domain
      vn = p0$indicators.variables[[iv]]
      sn = indicators.db( p=p0, DS="baseline", varnames=vn )
      yr_index = match( p$yrs, p0$yrs )
      yg = which(is.finite(yr_index))
      ym = which(is.na(yr_index))
      if (length(ym) > 0) {
        for ( vv in p$indicators.variables[[iv]] ) {
          PS[[vv]][yg] = sn[[vv]][,yr_index[yg]]
          PS[[vv]][ym] = rowMeans( sn[[vv]][], na.rm=TRUE )
        }
      } else {
        for ( vv in p$indicators.variables[[iv]] ) {
          PS[[vv]] = sn[[vv]][,yr_index]
        }
      }
    }

    ps_varnames = setdiff( p$varnames, p$variables$LOCS )

    PS = PS[ which(names(PS) %in% ps_varnames ) ] # time vars, if they are part of the model will be created within lbm

    oo = setdiff(p$varnames, c(ps_varnames, p$variables$LOCS) )
    if (length(oo) > 0 ) {
      print(oo )
      warning("Some variables are missing in the prediction surface, PS")
    }
    return (PS)
  }

  # -----------------------------------


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

        # qPS = quantile( PS, probs=p$lbm_quantile_bounds, na.rm=TRUE )
        # u = which( PS < qPS[1])
        # if (length(u)>0) PS[u] = qPS[1]
        # v = which( PS > qPS[2])
        # if (length(v)>0) PS[v] = qPS[2]
        
        # qPSsd = quantile( PSsd, probs=p$lbm_quantile_bounds, na.rm=TRUE )
        # u = which( PSsd < qPSsd[1])
        # if (length(u)>0) PSsd[u] = qPSsd[1]
        # v = which( PSsd > qPSsd[2])
        # if (length(v)>0) PSsd[v] = qPSsd[2]
      
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
        for (bvn in varnames ) {
          projectdir = file.path(p$project.root, "modelled", bvn, p$spatial.domain )
          outfile =  file.path( projectdir, paste( "snowcrab", "baseline", ret, p$spatial.domain, "rdata", sep= ".") )
          TS = NULL
          load( outfile)
          BL[[bvn]] = TS
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
        projectdir = file.path(p$project.root, "modelled", voi, p1$spatial.domain )
        dir.create( projectdir, recursive=T, showWarnings=F )

        L1 = bathymetry.db(p=p1, DS="baseline")
        nL1 = nrow(L1)
        
        TS = matrix( NA, nrow=nL1, ncol=p$ny )
        
        for (i in 1:p$ny ) {
          yr = p$yrs[i]
          TS[,i] = lbm_db( p=p1, DS="lbm.prediction", yr=yr, ret="mean")
         }

        outfile =  file.path( projectdir, paste( "snowcrab", "baseline", "mean", p1$spatial.domain, "rdata", sep= ".") )
        save( TS, file=outfile, compress=T )

        TS = matrix( NA, nrow=nL1, ncol=p$ny )
        for (i in 1:p$ny ) {
          yr = p$yrs[i]
          TS[,i] = lbm_db( p=p1, DS="lbm.prediction", yr=yr, ret="sd")
         }

        outfile =  file.path( projectdir, paste( "snowcrab", "baseline", "sd", p1$spatial.domain, "rdata", sep= ".") )
        save( TS, file=outfile, compress=T )
 
        print( outfile )

      }

      return( "Complete" )
 
    }


  # -----------------------

  if ( DS=="map.all" ) {

    allgrids = unique(c( p$spatial.domain.subareas, p$spatial.domain) )
    for ( gr in allgrids ) {
      print (gr)
      p1 = spatial_parameters(  p=p, type= gr )
      p1 = make.list( list( yrs=p1$yrs), Y=p1 )
      snowcrab_lbm( p=p1, DS="map.climatology" ) # no parallel option .. just a few
      parallel.run( snowcrab_lbm, p=p1, DS="map.annual", voi=voi )
    }

  }

  # -----------------------


  if ( DS %in% c("map.annual" ) ) {
    projectdir = file.path(p$project.root, "maps", voi, p$spatial.domain, "annual" )
    dir.create( projectdir, recursive=T, showWarnings=F )
    
    if (is.null(ip)) ip = 1:p$nruns
   
    # over-ride default dependent variable name if it exists
  
    loc = bathymetry.db(p=p, DS="baseline" )

    for (iy in ip ) {
      y = p$runs[iy, "yrs"]
      print(y)
      H = snowcrab_lbm( p=p, DS="predictions", year=y, ret="mean" )
      if (is.null(H)) next ()
      xyz = cbind(loc, H)
      uu = which( is.finite(rowSums(xyz)))
      if (length(uu) < 10) next()
      xyz = xyz[uu,]
      datarange = NULL
      datarange = snowcrab.lookup.mapparams( DS="datarange", voi ) # hardcoded data ranges 
      if (is.null(datarange)) {
        datarange=quantile(xyz[,3], probs=c(0.001,0.999), na.rm=TRUE) 
        datarange = seq( datarange[1], datarange[2], length.out=100 )
      }
      cols = color.code( "blue.black", datarange )
      annot = gsub( ".", " ", toupper(voi), fixed=TRUE )
      outfn = paste( voi, "mean", y, sep=".")
      lbm::lbm_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL, 
        loc=projectdir, fn=outfn, annot=annot, at=datarange , col.regions=cols,
        corners=p$corners, spatial.domain=p$spatial.domain ) 

      H = snowcrab_lbm( p=p, DS="predictions", year=y, ret="sd" )
      if (is.null(H)) next ()
      xyz = cbind(loc, H)
      uu = which( is.finite(rowSums(xyz)))
      if (length(uu) < 10) next()
      xyz = xyz[uu,]
      datarange = NULL      
      datarange = snowcrab.lookup.mapparams( DS="datarange", voi ) # hardcoded data ranges 
      if (is.null(datarange)) {
        datarange=quantile(xyz[,3], probs=c(0.001,0.999), na.rm=TRUE) 
        datarange = seq( datarange[1], datarange[2], length.out=100 )
      }
      cols = color.code( "blue.black", datarange )
      annot = gsub( ".", " ", toupper(voi), fixed=TRUE )
      outfn = paste( voi, "sd", y, sep=".")

      lbm::lbm_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL, 
        loc=projectdir, fn=outfn, annot=annot, at=datarange , col.regions=cols,
        corners=p$corners, spatial.domain=p$spatial.domain ) 
      print( file.path( projectdir, outfn))
    } 
    
  }


  # ------------------------------


  if ( DS %in% c("map.climatology" ) ) {
    projectdir = file.path(p$project.root, "maps", voi, p$spatial.domain, "climatology" )
    dir.create( projectdir, recursive=T, showWarnings=F )
    
    loc = bathymetry.db(p=p, DS="baseline" )

    H = snowcrab_lbm( p=p, DS="complete" )
    vnames = setdiff( names(H), c("plon", "plat" ))
    
    for (vn in vnames ) {
      xyz = cbind(loc, H[,vn])
      uu = which( is.finite(rowSums(xyz)))
      if (length(uu) < 10) next()
      xyz = xyz[uu,]
      datarange= NULL
      datarange = snowcrab.lookup.mapparams( DS="datarange", vn) # hardcoded data ranges 
      if (is.null(datarange)) {
        datarange=quantile(xyz[,3], probs=c(0.005,0.995), na.rm=TRUE) 
        datarange = seq( datarange[1], datarange[2], length.out=100 )
      }
      cols = color.code( "blue.black", datarange )
      annot = gsub( ".", " ", toupper(vn), fixed=TRUE )
      lbm::lbm_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL, 
        loc=projectdir, fn=vn, annot=annot, at=datarange, col.regions=cols,
        corners=p$corners, spatial.domain=p$spatial.domain ) 
      print( file.path( projectdir, vn))
    }

  }


}



