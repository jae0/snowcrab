
snowcrab_stmv = function( DS=NULL, p=NULL, year=NULL, ret="mean", varnames=NULL, ... ) {

  # deal with additional passed parameters
  # ---------------------
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ) ) 
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  # ---------------------
  # required:
  # over-ride default dependent variable name if it exists
  voi = NULL # variable of interest
  if (exists("selection",p)) if (exists("name", p$selection)) voi=p$selection$name
  if (is.null(voi)) if (exists("variables",p)) if (exists("Y", p$variables))    voi=p$variables$Y
  if (is.null(voi)) stop( "p$selection$name or p$variable$Y needs to be specified" )

  # ---------------------
  if (exists( "libs", p)) RLibrary( p$libs )
  if (!("stmv" %in% p$libs)) p$libs = c( p$libs, RLibrary( "stmv" ) )  # required for parallel processing

  # ---------------------
  if (DS=="parameters") {

    if (!exists("storage.backend", p)) p$storage.backend="bigmemory.ram"
    if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )

    if (!exists("boundary", p)) p$boundary = FALSE
    if (!exists("depth.filter", p)) p$depth.filter = 0 # depth (m) stats locations with elevation > 0 m as being on land (and so ignore)
    if (!exists("stmv_quantile_bounds", p)) p$stmv_quantile_bounds = c(0.025, 0.975) # remove these extremes in interpolations

    if (!exists("stmv_rsquared_threshold", p)) p$stmv_rsquared_threshold = 0.2 # lower threshold
    if (!exists("stmv_distance_statsgrid", p)) p$stmv_distance_statsgrid = 3 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    if (!exists("stmv_distance_prediction", p)) p$stmv_distance_prediction = p$stmv_distance_statsgrid*0.75  # this is a half window km
    if (!exists("stmv_distance_scale", p)) p$stmv_distance_scale = 50 # km ... approx guess of 95% AC range
    if (!exists("stmv_distance_min", p)) p$stmv_distance_min = 2
    if (!exists("stmv_distance_max", p)) p$stmv_distance_max = 65

    if (!exists("n.min", p)) p$n.min = 120 # n.min/n.max changes with resolution must be more than the number of knots/edf
    # min number of data points req before attempting to model timeseries in a localized space
    if (!exists("n.max", p)) p$n.max = 8000 # no real upper bound
    p$sampling = c( 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.25 )  #

    if (!exists("variables", p)) p$variables = list(
      Y = p$selection$name,
      LOCS = c("plon", "plat"),
      TIME = "tiyr",
      COV = c("z", "dZ", "ddZ", "log.substrate.grainsize", "t", "tmean.climatology", "tsd.climatology", "ca1", "ca2" ) )

    # additional variable to extract from aegis_db for inputs

    p$aegis_variables = list()
    # p$aegis_project_datasources = c("speciescomposition", "speciesarea", "sizespectrum", "condition", "metabolism", "biochem")
    p$aegis_project_datasources = "speciescomposition"
    for (id in p$aegis_project_datasources ) {
      pz = aegis::aegis_parameters( p=p, DS=id )
      pz_vars = intersect( pz$varstomodel, p$variables$COV )  # these are aegis vars to model
      if (length(pz_vars) > 0) p$aegis_variables[[id]] = pz_vars
    }

    if (!exists("stmv_variogram_method", p)) p$stmv_variogram_method = "fast"
    if (!exists("stmv_local_modelengine", p)) p$stmv_local_modelengine ="gam"
    if (!exists("stmv_global_modelengine", p)) p$stmv_global_modelengine ="gam"
    if (!exists("stmv_global_family", p)) p$stmv_global_family = gaussian(link=log)

    # using covariates as a first pass essentially makes it ~ kriging with external drift .. no time or space here
    if (!exists("stmv_global_modelformula", p)) p$stmv_global_modelformula = formula( paste(
      p$selection$name, ' ~ s(t, k=3, bs="ts") + s(tmean.climatology, k=3, bs="ts") + s(tsd.climatology, k=3, bs="ts")  ',
      ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
      ' + s(log.substrate.grainsize, k=3, bs="ts") + s(ca1, k=3, bs="ts") + s(ca2, k=3, bs="ts")   ' ))  # no space
    if (p$stmv_local_modelengine =="twostep") {
      # this is the time component (mostly) .. space enters as a rough constraint
      if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
        p$selection$name, '~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
          ' + s(cos.w, sin.w, yr, bs="ts", k=20) ',
          ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=20, bs="ts") ' ) )
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      # this is the spatial component
      # p$stmv_twostep_space = "spatial.process"
      # p$stmv_twostep_space = "fft"
      # p$stmv_twostep_space = "tps"
      if (!exists("stmv_twostep_space", p))  p$stmv_twostep_space = "krige"
      # if (!exists("stmv_twostep_space", p))  p$stmv_twostep_space = "tps"
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer=c("outer", "bfgs")
    }  else if (p$stmv_local_modelengine == "habitat") {
      p$stmv_global_family = binomial( link=log )
      if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
        p$selection$name, '~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
          ' + s(cos.w, sin.w, yr, bs="ts", k=10)  ',
          ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=10, bs="ts") ' ) )
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      # if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer="perf"
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer=c("outer", "bfgs")
    }  else if (p$stmv_local_modelengine == "gam") {
      if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
        p$selection$name, '~ s(yr, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
          ' + s(cos.w, sin.w, yr, bs="ts", k=25)  ',
          ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=25, bs="ts") ' ) )
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      # if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer="perf"
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer=c("outer", "bfgs")
    }  else if (p$stmv_local_modelengine == "bayesx") {
      # bayesx families are specified as characters, this forces it to pass as is and
      # then the next does the transformation internal to the "stmv__bayesx"
      # alternative models .. testing .. problem is that SE of fit is not accessible?
      p$stmv_local_modelformula = formula( paste(
        p$selection$name, ' ~ sx(yr, bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps") + sx(plon, bs="ps") + sx(plat,  bs="ps")',
          ' + sx(plon, plat, cos.w, sin.w, yr, bs="te") ' )
          # te is tensor spline
      )
      p$stmv_local_model_bayesxmethod="MCMC"
      p$stmv_local_model_distanceweighted = FALSE
    } else {
      message( "The specified stmv_local_modelengine is not tested/supported ... you are on your own ;) ..." )
    }
    return(p)
  }


  # --------------------------

  if (DS=="stmv_inputs") {
    # mostly based on aegis_db( DS="stmv_inputs")
    INP = snowcrab_stmv(p=p, DS="input_data" )  # , voi=p$selection$name
    PS  = snowcrab_stmv(p=p, DS="output_data"  ) # , voi=p$selection$name
    LOCS = bathymetry.db(p=p, DS="baseline")
    return (list(input=INP, output=list( LOCS=LOCS, COV=PS )) )
  }


  # --------------------------

  if (DS %in% c("input_data") ) {
    set = aegis::survey.db( p=p, DS="set.filter" ) # mature male > 95 mm

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

    coast = aegis::coastline.db( p=p, DS="mapdata.coastPolygon" )
    coast = spTransform( coast, CRS("+proj=longlat +datum=WGS84") )
    setcoord = SpatialPoints( as.matrix( set[, c("lon", "lat")]),  proj4string=CRS("+proj=longlat +datum=WGS84") )
    inside = sp::over( setcoord, coast )
    onland = which (is.finite(inside))
    if (length(onland)>0) set = set[-onland, ]

    set$lon = set$lat = NULL

    set$tiyr = lubridate::decimal_date( set$timestamp )

    bathy = bathymetry.db(p=p, DS="baseline")

    psse = spatial_parameters( p=p, spatial.domain = "SSE" )  # data are from this domain .. so far
    bathysse = bathymetry.db(p=psse, DS="baseline")

    # redo as set has changed
    locsmap = match(
      stmv::array_map( "xy->1", set[,c("plon","plat")], gridparams=p$gridparams ),
      stmv::array_map( "xy->1", bathy[,c("plon","plat")], gridparams=p$gridparams ) )

    # spatial vars and climatologies
    newvars = c("dZ", "ddZ", "log.substrate.grainsize", "tmean.climatology", "tsd.climatology", "b.range", "t.range" )
    sn = aegis_lookup( p=p, DS="spatial", locsmap=locsmap, varnames=newvars )
    set = cbind( set,  sn )

    oo = which( !is.finite(locsmap) )
    if (length(oo) > 0) {
      # catch stragglers from a larger domain
      locsmapsse = match(
        stmv::array_map( "xy->1", set[oo, c("plon","plat")], gridparams=psse$gridparams ),
        stmv::array_map( "xy->1", bathysse[,c("plon","plat")], gridparams=psse$gridparams ) )
      sn = aegis_lookup( p=psse, DS="spatial", locsmap=locsmapsse, varnames=newvars )
      for (nv in newvars) set[oo,nv] = sn[,nv]
    }

    # for space-time(year-averages)
    newvars = c( "tmean", "tsd", "tamplitude" ) 
    sn = aegis_lookup( p=p, DS="spatial.annual", locsmap=locsmap, timestamp=set[,"timestamp"], varnames=newvars )
    colnames( sn  ) = newvars
    set = cbind( set,  sn )

    nn = which( !is.finite(set$tmean) )
    if (length(nn) > 0) {
      # catch stragglers from a larger domain
      locsmapsse = match(
        stmv::array_map( "xy->1", set[nn, c("plon","plat")], gridparams=psse$gridparams ),
        stmv::array_map( "xy->1", bathysse[,c("plon","plat")], gridparams=psse$gridparams ) )
      sn = aegis_lookup( p=psse, DS="spatial.annual", locsmap=locsmapsse, timestamp=set[nn,"timestamp"], varnames=newvars )
      for (nv in newvars) set[nn,nv] = sn[,nv]
    }

    # names(set)[ names(set)=="amplitude"] ="tamplitude"

    # additional aegis_db variables
    for (iv in names(p$aegis_variables)) {
      p0 = aegis::aegis_parameters( p=p, DS=iv, year.assessment=p$year.assessment )
      p0 = aegis::spatial_parameters( p=p0, spatial.domain=p$spatial.domain ) # return to correct domain
      vn = p0$aegis_variables[[iv]]
      sn = aegis_lookup( p=p0, DS="spatial.annual", locsmap=locsmap, timestamp=set[,"timestamp"],
        varnames=vn, DB=aegis_db( p=p0, DS="baseline", varnames=vn ) )
      sn = as.data.frame(sn)
      names( sn  ) = p$aegis_variables[[iv]]
      set = cbind( set,  sn )

      mm = which( !is.finite(set[,vn[1]]) )
      if (length(mm) > 0) {
        # catch stragglers from a larger domain
        p0 = aegis::spatial_parameters( p=p0, spatial.domain="SSE" )
        locsmapsse = match(
          stmv::array_map( "xy->1", set[mm, c("plon","plat")], gridparams=p0$gridparams ),
          stmv::array_map( "xy->1", bathysse[,c("plon","plat")], gridparams=p0$gridparams ) )
        sn = aegis_lookup( p=p0, DS="spatial.annual", locsmap=locsmapsse, timestamp=set[mm,"timestamp"],
          varnames=vn, DB=aegis_db( p=p0, DS="baseline", varnames=vn ) )
        for (nv in vn) set[mm,nv] = sn[,nv]
      }
    }

    set = set[, which(names(set) %in% c( p$variables$LOCS, p$variables$COV, p$variables$Y, p$variables$TIME, "dyear", "yr",  "wt") ) ]  # a data frame
    oo = setdiff( c( p$variables$LOCS, p$variables$COV ), names(set))
    if (length(oo) > 0 ) {
      print(oo )
      warning("Some variables are missing in the input data")
    }
    set = na.omit(set)

    # cap quantiles of dependent vars
    dr = list()
    for (pvn in p$variables$COV) {
      dr[[pvn]] = quantile( set[,pvn], probs=p$stmv_quantile_bounds, na.rm=TRUE ) # use 95%CI
      il = which( set[,pvn] < dr[[pvn]][1] )
      if ( length(il) > 0 ) set[il,pvn] = dr[[pvn]][1]
      iu = which( set[,pvn] > dr[[pvn]][2] )
      if ( length(iu) > 0 ) set[iu,pvn] = dr[[pvn]][2]
    }

    return (set)

  }


  # ------------------------


  if (DS %in% c("output_data") ) {
    PS = aegis_db( p=p, DS="prediction.surface" ) # a list object with static and annually varying variables
    # names(PS)[ names(PS)=="amplitude"] ="tamplitude"

    # make years coherent for temperatures
    p0 = aegis::aegis_parameters(p=p, DS="temperature", year.assessment=p$year.assessment )
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

    # aegis_db variables
    for (iv in names(p$aegis_variables)) {
      p0 = aegis::aegis_parameters( p=p, DS=iv, year.assessment=p$year.assessment  )
      p0 = aegis::spatial_parameters( p=p0, spatial.domain=p$spatial.domain ) # return to correct domain

      vn = p0$aegis_variables[[iv]]
      sn = aegis_db( p=p0, DS="baseline", varnames=vn )
      yr_index = match( p$yrs, p0$yrs )
      yg = which(is.finite(yr_index))
      ym = which(is.na(yr_index))
      if (length(ym) > 0) {
        for ( vv in p$aegis_variables[[iv]] ) {
          PS[[vv]][yg] = sn[[vv]][,yr_index[yg]]
          PS[[vv]][ym] = rowMeans( sn[[vv]][], na.rm=TRUE )
        }
      } else {
        for ( vv in p$aegis_variables[[iv]] ) {
          PS[[vv]] = sn[[vv]][,yr_index]
        }
      }
    }

    PS = PS[ which(names(PS) %in% p$variables$COV ) ] # time vars, if they are part of the model will be created within stmv

    return (PS)
  }

  # -----------------------------------


  if ( DS %in% c("predictions", "predictions.redo" ) ) {
    # NOTE: the primary interpolated data were already created by stmv.
    # This routine points to this data and also creates
    # subsets of the data where required, determined by "spatial.domain.subareas"

    projectdir = file.path(p$data_root, "modelled", voi, p$spatial.domain )

    if (DS %in% c("predictions")) {
      P = V = W = NULL
      fn = file.path( projectdir, paste("stmv.prediction", ret,  year, "rdata", sep=".") )
      if (file.exists(fn) ) load(fn)
      if (ret=="mean") return (P)
      if (ret=="lb") return( V)
      if (ret=="ub") return( W)
    }

    parallel_run(
      p=p, 
      voi=voi, 
      runindex=list(yrs=p$yrs),
      FUNC = function( ip=NULL, p, voi ) {
        if (exists( "libs", p)) RLibrary( p$libs )
        if (is.null(ip)) ip = 1:p$nruns
        for ( ii in ip ) {
          # downscale and warp from p(0) -> p1
          year = p$runs[ii, "yrs"]
          # print (year)
          # default domain
          PP0 = stmv_db( p=p, DS="stmv.prediction", yr=year, ret="mean")
          VV0 = stmv_db( p=p, DS="stmv.prediction", yr=year, ret="lb")
          WW0 = stmv_db( p=p, DS="stmv.prediction", yr=year, ret="ub")
          p0 = spatial_parameters( p=p ) # from
          L0 = bathymetry.db( p=p0, DS="baseline" )
          L0i = stmv::array_map( "xy->2", L0[, c("plon", "plat")], gridparams=p0$gridparams )
          sreg = setdiff( p$spatial.domain.subareas, p$spatial.domain )

          for ( gr in sreg ) {
            p1 = spatial_parameters( p=p, spatial.domain=gr ) # 'warping' from p -> p1
            L1 = bathymetry.db( p=p1, DS="baseline" )
            L1i = stmv::array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
            L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
            L1$plon_1 = L1$plon # store original coords
            L1$plat_1 = L1$plat
            L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
            p1$wght = fields::setup.image.smooth( nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres, theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
            P = spatial_warp( PP0[], L0, L1, p0, p1, "fast", L0i, L1i )
            V = spatial_warp( VV0[], L0, L1, p0, p1, "fast", L0i, L1i )
            W = spatial_warp( WW0[], L0, L1, p0, p1, "fast", L0i, L1i )
            projectdir_p1 = file.path(p$data_root, "modelled", voi, p1$spatial.domain )
            dir.create( projectdir_p1, recursive=T, showWarnings=F )
            fn1_sg = file.path( projectdir_p1, paste("stmv.prediction.mean",  year, "rdata", sep=".") )
            fn2_sg = file.path( projectdir_p1, paste("stmv.prediction.lb",  year, "rdata", sep=".") )
            fn3_sg = file.path( projectdir_p1, paste("stmv.prediction.ub",  year, "rdata", sep=".") )
            save( P, file=fn1_sg, compress=T )
            save( V, file=fn2_sg, compress=T )
            save( W, file=fn3_sg, compress=T )
            print (fn1_sg)
          }
        }
        return()
      }
    )
    return ("Completed")

    if (0) {
      levelplot( P ~ plon_1 + plat_1, L1, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }

  }


  #  -------------------------------

  if (DS %in% c(  "stmv.stats", "stmv.stats.redo" )){


    if (DS %in% c("stmv.stats")) {
      stats = NULL
      projectdir = file.path(p$data_root, "modelled", voi, p$spatial.domain )
      fn = file.path( projectdir, paste( "stmv.statistics", "rdata", sep=".") )
      if (file.exists(fn) ) load(fn)
      return( stats )
    }

    sreg = setdiff( p$spatial.domain.subareas, p$spatial.domain )

    parallel_run(
      p=p, 
      voi=voi, 
      runindex=list(sreg=sreg),
      FUNC = function( ip=NULL, p, voi ) {
        if (exists( "libs", p)) RLibrary( p$libs )
        if (is.null(ip)) ip = 1:p$nruns
          # downscale and warp from p(0) -> p1 .. default domain
        S0 = stmv_db( p=p, DS="stats.to.prediction.grid" )
        Snames = colnames(S0) 
        p0 = spatial_parameters( p=p ) # from
        L0 = bathymetry.db( p=p0, DS="baseline" )
        L0i = stmv::array_map( "xy->2", L0[, c("plon", "plat")], gridparams=p0$gridparams )
        for ( ii in ip ) {
          gr = p$runs[ii,"sreg"]
          p1 = spatial_parameters( p=p, spatial.domain=gr ) # 'warping' from p -> p1
          L1 = bathymetry.db( p=p1, DS="baseline" )
          L1i = stmv::array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
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
          projectdir_p1 = file.path(p$data_root, "modelled", voi, p1$spatial.domain )
          dir.create( projectdir_p1, recursive=T, showWarnings=F )
          fn1_sg = file.path( projectdir_p1, paste("stmv.statistics", "rdata", sep=".") )
          save( stats, file=fn1_sg, compress=T )
          print (fn1_sg)
        }
      }
    )
  
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
      projectdir = file.path(p$data_root, "modelled", voi, p$spatial.domain )
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

      p1 = spatial_parameters( p=p, spatial.domain=gr ) #target projection
      L1 = bathymetry.db(p=p1, DS="baseline")

      BS = snowcrab_stmv( p=p1, DS="stmv.stats" )
      colnames(BS) = paste(voi, colnames(BS), sep=".")
      IC = cbind( L1, BS )

      # climatology
      nL1 = nrow(L1)
      PS = PSlb = PSub = matrix( NA, nrow=nL1, ncol=p$ny )
      if (is.null(voi)) p1$variables$Y = voi # need to send this to get the correct results
      for (iy in 1:p$ny) {
        yr = p$yrs[iy]
        PS[,iy] = stmv_db( p=p1, DS="stmv.prediction", yr=yr, ret="mean")
        PSlb[,iy] = stmv_db( p=p1, DS="stmv.prediction", yr=yr, ret="lb")
        PSub[,iy] = stmv_db( p=p1, DS="stmv.prediction", yr=yr, ret="ub")
      }

      # qPS = quantile( PS, probs=p$stmv_quantile_bounds, na.rm=TRUE )
      # u = which( PS < qPS[1])
      # if (length(u)>0) PS[u] = qPS[1]
      # v = which( PS > qPS[2])
      # if (length(v)>0) PS[v] = qPS[2]

      # qPSlb = quantile( PSlb, probs=p$stmv_quantile_bounds, na.rm=TRUE )
      # u = which( PSlb < qPSlb[1])
      # if (length(u)>0) PSlb[u] = qPSlb[1]
      # v = which( PSlb > qPSlb[2])
      # if (length(v)>0) PSlb[v] = qPSlb[2]

      # qPSub = quantile( PSub, probs=p$stmv_quantile_bounds, na.rm=TRUE )
      # u = which( PSub < qPSub[1])
      # if (length(u)>0) PSub[u] = qPSub[1]
      # v = which( PSub > qPSub[2])
      # if (length(v)>0) PSub[v] = qPSub[2]

      CL = cbind( apply( PS, 1, mean, na.rm=TRUE ),
                  apply( PSlb, 1, mean, na.rm=TRUE ),
                  apply( PSub, 1, mean, na.rm=TRUE ) )
      colnames(CL) = paste( voi, c("mean", "lb", "ub"), "climatology", sep=".")
      IC = cbind( IC, CL )
      PS = PSlb = PSub = NULL

      projectdir = file.path(p$data_root, "modelled", voi, p1$spatial.domain )
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
        projectdir = file.path(p$data_root, "modelled", bvn, p$spatial.domain )
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


    parallel_run(
      p=p, 
      voi=voi, 
      runindex=list(grids=grids),
      FUNC = function( ip=NULL, p, voi ) {
        if (exists( "libs", p)) RLibrary( p$libs )
        if (is.null(ip)) ip = 1:p$nruns
        
        for (ii in ip ) {
          gr = p$runs[ii,"grids"]
          print(gr)
          p1 = spatial_parameters( p=p, spatial.domain=gr ) #target projection
          projectdir = file.path(p$data_root, "modelled", voi, p1$spatial.domain )
          dir.create( projectdir, recursive=T, showWarnings=F )
          L1 = bathymetry.db(p=p1, DS="baseline")
          nL1 = nrow(L1)
          TS = matrix( NA, nrow=nL1, ncol=p$ny )
          for (i in 1:p$ny ) {
            TS[,i] = stmv_db( p=p1, DS="stmv.prediction", yr=p$yrs[i], ret="mean")
          }
          outfile =  file.path( projectdir, paste( "snowcrab", "baseline", "mean", p1$spatial.domain, "rdata", sep= ".") )
          save( TS, file=outfile, compress=T )
          TS = TS[] * NA
          for (i in 1:p$ny ) {
            TS[,i] = stmv_db( p=p1, DS="stmv.prediction", yr=p$yrs[i], ret="lb")
          }
          outfile =  file.path( projectdir, paste( "snowcrab", "baseline", "lb", p1$spatial.domain, "rdata", sep= ".") )
          save( TS, file=outfile, compress=T )
          TS = TS[] * NA
          for (i in 1:p$ny ) {
            TS[,i] = stmv_db( p=p1, DS="stmv.prediction", yr=p$yrs[i], ret="ub")
          }
          outfile =  file.path( projectdir, paste( "snowcrab", "baseline", "ub", p1$spatial.domain, "rdata", sep= ".") )
          save( TS, file=outfile, compress=T )
        }
        return(NULL)
      }
    )
    return( "Complete" )
  }


  # -----------------------

  if ( DS=="map.all" ) {

    allgrids = unique(c( p$spatial.domain.subareas, p$spatial.domain) )
    for ( gr in allgrids ) {
      print (gr)
      p1 = spatial_parameters(  p=p, spatial.domain= gr )
      snowcrab_stmv( p=p1, DS="map.climatology", voi=voi ) # no parallel option .. just a few
      snowcrab_stmv( p=p1, DS="map.annual", voi=voi ) 
    }

  }

  # -----------------------


  if ( DS %in% c("map.annual" ) ) {

    if (exists( "libs", p)) RLibrary( p$libs )

    for ( year in p$yrs ) {
        projectdir = file.path(p$data_root, "maps", voi, p$spatial.domain, "annual" )
        dir.create( projectdir, recursive=T, showWarnings=F )
        loc = bathymetry.db(p=p, DS="baseline" )

        # downscale and warp from p(0) -> p1
        # print(year)
        H = snowcrab_stmv( p=p, DS="predictions", year=year, ret="mean" )
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
        outfn = paste( voi, "mean", year, sep=".")
        
        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL,
          annot=annot, at=datarange, col.regions=cols,
          corners=p$corners, spatial.domain=p$spatial.domain )
        print(lp)
        dev.off()
        
        H = snowcrab_stmv( p=p, DS="predictions", year=year, ret="lb" )
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
        outfn = paste( voi, "lb", year, sep=".")
        
        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL,
          annot=annot, at=datarange, col.regions=cols,
          corners=p$corners, spatial.domain=p$spatial.domain )
        print(lp)
        dev.off()

        H = snowcrab_stmv( p=p, DS="predictions", year=year, ret="ub" )
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
        outfn = paste( voi, "ub", year, sep=".")
        
        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL,
              annot=annot, at=datarange , col.regions=cols,
              corners=p$corners, spatial.domain=p$spatial.domain )
        print(lp)
        dev.off()

        print( file.path( projectdir, outfn))
    }
    return("Finished")
  }


  # ------------------------------


  if ( DS %in% c("map.climatology" ) ) {
    if (exists( "libs", p)) RLibrary( p$libs )
    H = snowcrab_stmv( p=p, DS="complete" )
    vnames = setdiff( names(H), c("plon", "plat" ))
    H = NULL
 
    for (vn in vnames) {
        projectdir = file.path(p$data_root, "maps", voi, p$spatial.domain, "climatology" )
        dir.create( projectdir, recursive=T, showWarnings=F )
        loc = bathymetry.db(p=p, DS="baseline" )
        H = snowcrab_stmv( p=p, DS="complete" )
        vnames = setdiff( names(H), c("plon", "plat" ))

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

        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(vn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL,
          annot=annot, at=datarange, col.regions=cols,
          corners=p$corners, spatial.domain=p$spatial.domain )
        print(lp)
        dev.off()

        print( file.path( projectdir, vn))
    }
    return( "Completed")
  }
}
