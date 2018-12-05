
snowcrab_stmv = function( DS=NULL, p=NULL, year=NULL, ret="mean", varnames=NULL, ... ) {

  # deal with additional passed parameters
  # ---------------------
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ) )
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  # ---------------------
  if (exists( "libs", p)) RLibrary( p$libs )
  if (!("stmv" %in% p$libs)) p$libs = c( p$libs, RLibrary( "stmv" ) )  # required for parallel processing

  # due to formulae being created on the fly, these are required params
  # if (!exists("variables", p)) stop("Please define p$variables$Y")
  # if (!exists("Y", p$variables)) stop("Please define p$variables$Y")

  if (!exists("variables", p)) p$variables = list()

  if (!exists("LOCS", p$variables)) p$variables$LOCS = c("plon", "plat")
  if (!exists("TIME", p$variables)) p$variables$TIME = "tiyr"

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

    if (!exists("n.min", p)) p$n.min = 100 # n.min/n.max changes with resolution must be more than the number of knots/edf
    # min number of data points req before attempting to model timeseries in a localized space
    if (!exists("n.max", p)) p$n.max = 6000 # actually can have a lot of data from logbooks ... this keeps things reasonable in terms of run-time
    p$sampling = c( 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.25 )  #

    p$stmv_dimensionality="space-year"

    # due to formulae being potentially created on the fly, these are required params
    if (!exists("Y", p$variables)) {
      if (exists("stmv_local_modelformula", p))  {
        if (!is.null(p$stmv_local_modelformula)) {
          if (p$stmv_local_modelformula != "none") {
            oo = all.vars( p$stmv_local_modelformula[[2]] )
            if (length(oo) > 0) p$variables$Y = oo
          }
        }
      }
      if (exists("stmv_global_modelformula", p))  {
        if (!is.null(p$stmv_global_modelformula)) {
          if (p$stmv_global_modelformula != "none") {
            oo = all.vars( p$stmv_global_modelformula[[2]] )
            if (length(oo) > 0) p$variables$Y = oo
          }
        }
      }
    }
    if (!exists("Y", p$variables)) p$variables$Y = "not_defined" # this can be called to get covars.. do not stop


    # additional variable to extract from aegis_db for inputs
    p$aegis_variables = list()
    # p$aegis_project_datasources = c("speciescomposition", "speciesarea", "sizespectrum", "condition", "metabolism", "biochem")
    if (!exists("aegis_project_datasources", p)) p$aegis_project_datasources = "speciescomposition"
    for (id in p$aegis_project_datasources ) {
      pz = aegis::aegis_parameters( p=p, DS=id )
      pz_vars = intersect( pz$varstomodel, p$variables$COV )  # these are aegis vars to model
      if (length(pz_vars) > 0) p$aegis_variables[[id]] = pz_vars
    }

    if (!exists("stmv_variogram_method", p)) p$stmv_variogram_method = "gstat"
    if (!exists("stmv_local_modelengine", p)) p$stmv_local_modelengine ="gam"

    if (!exists("stmv_global_modelengine", p)) p$stmv_global_modelengine ="gam"
    if (!exists("stmv_global_family", p)) p$stmv_global_family = gaussian(link="log")

    # using covariates as a first pass essentially makes it ~ kriging with external drift .. no time or space here
    if (!exists("stmv_global_modelformula", p)) p$stmv_global_modelformula = formula( paste(
      p$variables$Y, ' ~ s(t, k=3, bs="ts") + s(tmean.climatology, k=3, bs="ts") + s(tsd.climatology, k=3, bs="ts")  ',
      ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
      ' + s(log(substrate.grainsize), k=3, bs="ts") + s(pca1, k=3, bs="ts") + s(pca2, k=3, bs="ts")   ' ))  # no space

    if (p$stmv_local_modelengine =="twostep") {
      # this is the time component (mostly) .. space enters as a rough constraint
      # if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
      #   p$variables$Y, '~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
      #     ' + s(cos.w, sin.w, yr, bs="ts", k=20) ',
      #     ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=20, bs="ts") ' ) )
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer=c("outer", "bfgs")

      if (!exists("stmv_twostep_time", p))  p$stmv_twostep_time = "gam"
      if (p$stmv_twostep_time == "gam") {
        if (!exists("stmv_local_modelformula_time", p))  p$stmv_local_modelformula_time = formula( paste(
          p$variables$Y, '~ s(yr, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
            ' + s(cos.w, sin.w, yr, bs="ts", k=25)  ',
            ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=25, bs="ts") ' ) )
      }

      # this is the spatial component
      # p$stmv_twostep_space = "spatial.process"
      # p$stmv_twostep_space = "fft"
      # p$stmv_twostep_space = "tps"
      # if (!exists("stmv_twostep_space", p))  p$stmv_twostep_space = "krige"
      # if (!exists("stmv_twostep_space", p))  p$stmv_twostep_space = "tps"

      if (!exists("stmv_twostep_space", p))  p$stmv_twostep_space = "fft" #  fft==spatial.process, krige (very slow), lowpass, lowpass_spatial.process
      if (p$stmv_twostep_space == "gam") {
        if (!exists("stmv_local_modelformula_space", p))  p$stmv_local_modelformula_space = formula( paste(
        p$variables$Y, '~ s(log(z), k=3, bs="ts") + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s( log(z), plon, plat, k=27, bs="ts")  ') )
      }
      if (!exists("stmv_fft_filter", p)) p$stmv_fft_filter="spatial.process" #  fft==spatial.process, krige (very slow), lowpass, lowpass_spatial.process


    }  else if (p$stmv_local_modelengine == "gam") {

      if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
        p$variables$Y, '~ s(yr, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
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
        p$variables$Y, ' ~ sx(yr, bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps") + sx(plon, bs="ps") + sx(plat,  bs="ps")',
          ' + sx(plon, plat, cos.w, sin.w, yr, bs="te") ' )
          # te is tensor spline
      )
      p$stmv_local_model_bayesxmethod="MCMC"
      p$stmv_local_model_distanceweighted = FALSE
    } else {
      message( "The specified stmv_local_modelengine is not tested/supported ... you are on your own ;) ..." )
    }

    p = stmv_variablelist(p=p)  # decompose into covariates, etc

    return(p)
  }


  # --------------------------

  if (DS=="stmv_inputs") {
    # mostly based on aegis_db( DS="stmv_inputs")
    INP = snowcrab_stmv(p=p, DS="input_data" )
    PS  = snowcrab_stmv(p=p, DS="output_data" )
    # alternatively using aegis:
    # PS = aegis_db_extract( vars=p$variables$COV, yrs=p$yrs, spatial.domain=p$spatial.domain, dyear=p$prediction.dyear )
    LOCS = bathymetry.db(p=p, DS="baseline")
    return (list(input=INP, output=list( LOCS=LOCS, COV=PS )) )
  }


  # --------------------------

  if (DS %in% c("input_data") ) {
    set = aegis::survey.db( p=p, DS="det.filter" ) # mature male > 95 mm

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
      names(set)[ which( names(set) =="totmass")] = p$variables$Y
      set$Y = NULL
      set$wt = set$sa
    }

    if ( p$selection$type=="presence_absence") {
      # must run here as we need the wgt from this for both PA and abundance
      if (p$selection$survey$drop.groundfish.data) {
          todrop = which(set$data.source == "groundfish")
          if (length(todrop)> 0 ) set = set[ -todrop, ]
      }

      if ( grepl( "snowcrab.large.males", p$variables$Y ) ) {
        # add commerical fishery data --
        # depth data is problematic ... drop for now
        lgbk = logbook.db( DS="fisheries.complete", p=p )
        lgbk = lgbk[ which( is.finite( lgbk$landings)), ]
        lgbk = lgbk[ which( lgbk$year > 2005), ]  # previous to this all sorts of traps were used
        lgbk = lgbk[ which( as.numeric(lgbk$soak.time) >= 12 & as.numeric(lgbk$soak.time) <= 48), ]   # avoid nonlinearity in catch with time
        lgbk$zm = lgbk$cpue / as.numeric(lgbk$soak.time)  # approx with realtive catch rate in time
        lgbk$totmass = NA # dummy to bring in mass as well
        lgbk$data.source = "logbooks"
        lgbk$z = exp( lgbk$z )
        nms = intersect( names(set) , names( lgbk) )
        set = rbind( set[, nms], lgbk[,nms] )
      }

      set = presence.absence( X=set, vname="zm", px=p$habitat.threshold.quantile )  # determine presence absence and weighting
      names(set)[ which( names(set) =="Y")] = p$variables$Y
      set$totmass = NULL
      set = set[ which(is.finite(set$plon + set$plat)),]
    }

    set = set[ which(is.finite(set[, p$variables$Y])),]

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

    # redo as "set" has changed
    locsmap = match(
      stmv::array_map( "xy->1", set[,c("plon","plat")], gridparams=p$gridparams ),
      stmv::array_map( "xy->1", bathy[,c("plon","plat")], gridparams=p$gridparams ) )

    # spatial vars and climatologies
    # newvars = c("dZ", "ddZ", "substrate.grainsize", "tmean.climatology", "tsd.climatology", "b.range", "t.range" )
    newvars = setdiff(p$variables$COV, names(set) )
    if (length(newvars) > 0) {
      sn = NULL
      sn = aegis_lookup( p=p, DS="spatial", locsmap=locsmap, varnames=newvars )
      if (!is.null(sn)) {
        set = cbind( set,  sn )
      }
    }

    oo = which( !is.finite(locsmap) )
    if (length(oo) > 0) {
      # catch stragglers from a larger domain
      locsmapsse = match(
        stmv::array_map( "xy->1", set[oo, c("plon","plat")], gridparams=psse$gridparams ),
        stmv::array_map( "xy->1", bathysse[,c("plon","plat")], gridparams=psse$gridparams ) )
        if (length(newvars) > 0) {
          sn = NULL
          sn = aegis_lookup( p=psse, DS="spatial", locsmap=locsmapsse, varnames=newvars )
          if (!is.null(sn)) {
            for (nv in names(sn)) {
              if (is.vector(sn)) {
                set[oo,nv] = sn[nv]
              } else {
                set[oo,nv] = sn[,nv]
              }
            }
          }
        }
    }

    # for space-time(year-averages)
    # newvars = c( "tmean", "tsd", "tamplitude" )
    newvars = setdiff(p$variables$COV, names(set) )
    if (length(newvars) > 0) {
      sn = NULL
      sn = aegis_lookup( p=p, DS="spatial.annual", locsmap=locsmap, timestamp=set[,"timestamp"], varnames=newvars )
      if (!is.null(sn)) {
        colnames( sn  ) = newvars
        set = cbind( set,  sn )
      }
    }

    nn = which( !is.finite(set$tmean) )
    if (length(nn) > 0) {
      # catch stragglers from a larger domain
      locsmapsse = match(
        stmv::array_map( "xy->1", set[nn, c("plon","plat")], gridparams=psse$gridparams ),
        stmv::array_map( "xy->1", bathysse[,c("plon","plat")], gridparams=psse$gridparams ) )
      sn = aegis_lookup( p=psse, DS="spatial.annual", locsmap=locsmapsse, timestamp=set[nn,"timestamp"], varnames=newvars )
      for (nv in newvars) set[nn,nv] = sn[,nv]
    }



    # for space-time(year-averages)
    # newvars = c( "tmean", "tsd", "tamplitude" )
    newvars = setdiff(p$variables$COV, names(set) )
    if (length(newvars) > 0) {
      sn = NULL
      sn = aegis_lookup( p=p, DS="spatial.annual.seasonal", locsmap=locsmap, timestamp=set[,"timestamp"], varnames=newvars )
      if (!is.null(sn)) {
        colnames( sn  ) = newvars
        set = cbind( set,  sn )
      }
    }

    nn = which( !is.finite(set$tmean) )
    if (length(nn) > 0) {
      # catch stragglers from a larger domain
      locsmapsse = match(
        stmv::array_map( "xy->1", set[nn, c("plon","plat")], gridparams=psse$gridparams ),
        stmv::array_map( "xy->1", bathysse[,c("plon","plat")], gridparams=psse$gridparams ) )
      sn = aegis_lookup( p=psse, DS="spatial.annual.seasonal", locsmap=locsmapsse, timestamp=set[nn,"timestamp"], varnames=newvars )
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

    # alternative to below is... still in testing
    # PS = aegis_db_extract( varnames=p$variables$COV, yrs=p$yrs, spatial.domain=p$spatial.domain, aegis_project_datasources =p$aegis_variables )

    PS = aegis_db( p=p, DS="prediction.surface" ) # a list object with static and annually varying variables
    # names(PS)[ names(PS)=="amplitude"] ="tamplitude"

    # make years coherent for temperatures
    PSyrs = colnames(PS[["t"]])
    pt = aegis::aegis_parameters(p=p, DS="temperature")
    yr_index = match( as.character(p$yrs), PSyrs )
    yg = which(is.finite(yr_index))
    ym = which(is.na(yr_index))
    if (length(ym) > 0) {
      for ( vn in c("t", pt$bstats) ) {  # annual summaries
        PS[[vn]][yg] = PS[[vn]][,yr_index[yg]]
        PS[[vn]][ym] = rowMeans( PS[[vn]][], na.rm=TRUE )
      }
    } else {
      for ( vn in c("t", pt$bstats) ) {
        PS[[vn]][] = PS[[vn]][,yr_index]
      }
    }

    # aegis_db variables
    for (iv in names(p$aegis_variables)) {
      pv = aegis::aegis_parameters( p=p, DS=iv, year.assessment=p$year.assessment  )
      pv = aegis::spatial_parameters( p=pv, spatial.domain=p$spatial.domain ) # return to correct domain

      vn = pv$aegis_variables[[iv]]
      sn = aegis_db( p=pv, DS="baseline", varnames=vn )
      snyrs = colnames(sn[[vn[1]]])

      yr_index = match( as.character(p$yrs), snyrs )
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
    # not strictly required for snow crab data analysis as there are no sub-areas that are analysed
    # at present, but in case there are small-area analysis in future, this is a mechnanism

    projectdir = file.path(p$data_root, "modelled", p$variables$Y, p$spatial.domain )

    if (DS %in% c("predictions")) {
      P = Pl = Pu = NULL
      fn = file.path( projectdir, paste("stmv.prediction", ret,  year, "rdata", sep=".") )
      if (file.exists(fn) ) load(fn)
      if (ret=="mean") return (P)
      if (ret=="lb") return( Pl)
      if (ret=="ub") return( Pu)
    }

    sreg = setdiff( p$spatial.domain.subareas, p$spatial.domain ) # see  note above
    if (is.null(sreg)) return
    if (length(sreg) < 1 ) return

    p0 = spatial_parameters( p=p ) # make explicit
    L0 = bathymetry.db( p=p0, DS="baseline" )
    L0i = stmv::array_map( "xy->2", L0[, c("plon", "plat")], gridparams=p0$gridparams )

    for ( year in p$yrs ) {
      # print (year)
      # default domain
      PP0 = stmv_db( p=p, DS="stmv.prediction", yr=year, ret="mean")
      VV0 = stmv_db( p=p, DS="stmv.prediction", yr=year, ret="lb")
      WW0 = stmv_db( p=p, DS="stmv.prediction", yr=year, ret="ub")

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
        Pl = spatial_warp( VV0[], L0, L1, p0, p1, "fast", L0i, L1i )
        Pu = spatial_warp( WW0[], L0, L1, p0, p1, "fast", L0i, L1i )
        projectdir_p1 = file.path(p$data_root, "modelled", p1$variables$Y, p1$spatial.domain )
        dir.create( projectdir_p1, recursive=T, showWarnings=F )
        fn1_sg = file.path( projectdir_p1, paste("stmv.prediction.mean",  year, "rdata", sep=".") )
        fn2_sg = file.path( projectdir_p1, paste("stmv.prediction.lb",  year, "rdata", sep=".") )
        fn3_sg = file.path( projectdir_p1, paste("stmv.prediction.ub",  year, "rdata", sep=".") )
        save( P, file=fn1_sg, compress=T )
        save( Pl, file=fn2_sg, compress=T )
        save( Pu, file=fn3_sg, compress=T )
        print (fn1_sg)
      }
    }

    return ("Completed")

    if (0) {
      levelplot( P ~ plon_1 + plat_1, L1, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }

  }


  #  -------------------------------

  if (DS %in% c(  "stmv.stats", "stmv.stats.redo" )){

    if (DS %in% c("stmv.stats")) {
      stats = NULL
      projectdir = file.path(p$data_root, "modelled", p$variables$Y, p$spatial.domain )
      fn = file.path( projectdir, paste( "stmv.statistics", "rdata", sep=".") )
      if (file.exists(fn) ) load(fn)
      return( stats )
    }

    sreg = setdiff( p$spatial.domain.subareas, p$spatial.domain )
    if (is.null(sreg)) return
    if (length(sreg) < 1 ) return

    S0 = stmv_db( p=p, DS="stmv.stats" )
    Snames = colnames(S0)
    p0 = spatial_parameters( p=p ) # from
    L0 = bathymetry.db( p=p0, DS="baseline" )
    L0i = stmv::array_map( "xy->2", L0[, c("plon", "plat")], gridparams=p0$gridparams )

    for ( gr in sreg ) {
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
      projectdir_p1 = file.path(p$data_root, "modelled", p$variables$Y, p1$spatial.domain )
      dir.create( projectdir_p1, recursive=T, showWarnings=F )
      fn1_sg = file.path( projectdir_p1, paste("stmv.statistics", "rdata", sep=".") )
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
    # assemble data for use by other projects
    if (DS=="complete") {
      IC = NULL
      projectdir = file.path(p$data_root, "modelled", p$variables$Y, p$spatial.domain )
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

    grids = unique( c(p$spatial.domain.subareas , p$spatial.domain ) ) # operate upon every domain

    # if ( p$variables$Y=="snowcrab.large.males_abundance" ) {
    #   # copied from  snowcrab_stmv(p=p, DS="input_data" )
    #   set = set[ which(set$data.source == "snowcrab"), ]
    #   qq = quantile( set$totmass, probs=0.975, na.rm=TRUE )
    #   qb = c( 0, qq )   # emprical range limit (do not extrapolate)
    # }

    for (gr in grids ) {
      p1 = spatial_parameters( p=p, spatial.domain=gr ) #target projection
      # p1$variables$Y = p$variables$Y # need to send this to get the correct results
      L1 = bathymetry.db(p=p1, DS="baseline")
      BS = snowcrab_stmv( p=p1, DS="stmv.stats" )
      colnames(BS) = paste(p$variables$Y, colnames(BS), sep=".")
      IC = cbind( L1, BS )
      # climatology
      nL1 = nrow(L1)
      PS = PSlb = PSub = matrix( NA, nrow=nL1, ncol=p$ny )
      for (iy in 1:p$ny) {
        yr = p$yrs[iy]
        PS[,iy] = stmv_db( p=p1, DS="stmv.prediction", yr=yr, ret="mean")
        PSlb[,iy] = stmv_db( p=p1, DS="stmv.prediction", yr=yr, ret="lb")
        PSub[,iy] = stmv_db( p=p1, DS="stmv.prediction", yr=yr, ret="ub")
      }
      CL = cbind( apply( PS, 1, mean, na.rm=TRUE ),
                  apply( PSlb, 1, mean, na.rm=TRUE ),
                  apply( PSub, 1, mean, na.rm=TRUE ) )
      colnames(CL) = paste( p1$variables$Y, c("mean", "lb", "ub"), "climatology", sep=".")
      IC = cbind( IC, CL )
      PS = PSlb = PSub = NULL
      projectdir = file.path(p$data_root, "modelled", p1$variables$Y, p1$spatial.domain )
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

    grids = unique( c(p$spatial.domain.subareas , p$spatial.domain ) ) # operate upon every domain

    for (gr in grids) {
        print(gr)
        p1 = spatial_parameters( p=p, spatial.domain=gr ) #target projection
        projectdir = file.path(p$data_root, "modelled", p$variables$Y, p1$spatial.domain )
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
    return( "Complete" )
  }


  # -----------------------

  if ( DS=="map.all" ) {

    allgrids = unique(c( p$spatial.domain.subareas, p$spatial.domain) )
    for ( gr in allgrids ) {
      print (gr)
      p1 = spatial_parameters(  p=p, spatial.domain= gr )
      snowcrab_stmv( p=p1, DS="map.climatology" ) # no parallel option .. just a few
      snowcrab_stmv( p=p1, DS="map.annual" )
    }

  }

  # -----------------------


  if ( DS %in% c("map.annual" ) ) {

    annot.cex=0.65
    eps = 0.001

    for ( year in p$yrs ) {
        projectdir = file.path(p$data_root, "maps", p$variables$Y, p$spatial.domain, "annual" )
        dir.create( projectdir, recursive=T, showWarnings=F )
        loc = bathymetry.db(p=p, DS="baseline" )

        # downscale and warp from p(0) -> p1
        # print(year)
        H = snowcrab_stmv( p=p, DS="predictions", year=year, ret="mean" )
        if (is.null(H)) next ()
        # H = log(H)
        xyz = cbind(loc, H)
        uu = which( is.finite(rowSums(xyz)))
        if (length(uu) < 10) next()
        xyz = xyz[uu,]
        datarange = NULL
        datarange = snowcrab.lookup.mapparams( DS="datarange", p$variables$Y ) # hardcoded data ranges
        if (is.null(datarange)) {
          datarange=quantile(xyz[,3], probs=c(0.001,0.999), na.rm=TRUE)
          datarange = seq( datarange[1], datarange[2], length.out=100 )
        }
        cols = color.code( "blue.black", datarange )
        annot = gsub( ".", " ", toupper(p$variables$Y), fixed=TRUE )
        outfn = paste( p$variables$Y, "mean", year, sep=".")

        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL,
          annot=annot, annot.cex=annot.cex, at=datarange, col.regions=cols,
          corners=p$corners, spatial.domain=p$spatial.domain )
        print(lp)
        dev.off()

        H = snowcrab_stmv( p=p, DS="predictions", year=year, ret="lb" )
        if (is.null(H)) next ()
        # H = log(H)
        xyz = cbind(loc, H)
        uu = which( is.finite(rowSums(xyz)))
        if (length(uu) < 10) next()
        xyz = xyz[uu,]
        datarange = NULL
        datarange = snowcrab.lookup.mapparams( DS="datarange", p$variables$Y ) # hardcoded data ranges
        if (is.null(datarange)) {
          datarange=quantile(xyz[,3], probs=c(0.001,0.999), na.rm=TRUE)
          if (diff(datarange) < eps) datarange[2] = datarange[2]+ eps
          datarange = seq( datarange[1], datarange[2], length.out=100 )
        }
        cols = color.code( "blue.black", datarange )
        annot = gsub( ".", " ", toupper(p$variables$Y), fixed=TRUE )
        outfn = paste( p$variables$Y, "lb", year, sep=".")

        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL,
          annot=annot, annot.cex=annot.cex, at=datarange, col.regions=cols,
          corners=p$corners, spatial.domain=p$spatial.domain )
        print(lp)
        dev.off()

        H = snowcrab_stmv( p=p, DS="predictions", year=year, ret="ub" )
        if (is.null(H)) next ()
        # H = log(H)
        xyz = cbind(loc, H)
        uu = which( is.finite(rowSums(xyz)))
        if (length(uu) < 10) next()
        xyz = xyz[uu,]
        datarange = NULL
        datarange = snowcrab.lookup.mapparams( DS="datarange", p$variables$Y ) # hardcoded data ranges
        if (is.null(datarange)) {
          datarange=quantile(xyz[,3], probs=c(0.001,0.999), na.rm=TRUE)
          datarange = seq( datarange[1], datarange[2], length.out=100 )
        }
        cols = color.code( "blue.black", datarange )
        annot = gsub( ".", " ", toupper(p$variables$Y), fixed=TRUE )
        outfn = paste( p$variables$Y, "ub", year, sep=".")

        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL,
              annot=annot, annot.cex=annot.cex, at=datarange , col.regions=cols,
              corners=p$corners, spatial.domain=p$spatial.domain )
        print(lp)
        dev.off()

        print( file.path( projectdir, outfn))
    }
    return("Finished")
  }


  # ------------------------------


  if ( DS %in% c("map.climatology" ) ) {

    annot.cex=0.75
    eps = 0.001

    H = snowcrab_stmv( p=p, DS="complete" )
    vnames = setdiff( names(H), c("plon", "plat" ))
    H = NULL

    for (vn in vnames) {
        projectdir = file.path(p$data_root, "maps", p$variables$Y, p$spatial.domain, "climatology" )
        dir.create( projectdir, recursive=T, showWarnings=F )
        loc = bathymetry.db(p=p, DS="baseline" )
        H = snowcrab_stmv( p=p, DS="complete" )
        vnames = setdiff( names(H), c("plon", "plat" ))

        xyz = cbind(loc, H[,vn])
        # if (grepl("abundance", vn)) xyz[,3] = log(xyz[,3])
        uu = which( is.finite(rowSums(xyz)))
        if (length(uu) < 10) next()
        xyz = xyz[uu,]
        datarange= NULL
        datarange = snowcrab.lookup.mapparams( DS="datarange", vn) # hardcoded data ranges
        if (is.null(datarange)) {
          datarange=quantile(xyz[,3], probs=c(0.005,0.995), na.rm=TRUE)
          if (diff(datarange) < eps) datarange[2] = datarange[2]+ eps
          datarange = seq( datarange[1], datarange[2], length.out=100 )
        }
        cols = color.code( "blue.black", datarange )
        annot = gsub( ".", " ", toupper(vn), fixed=TRUE )

        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(vn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=FALSE, depthcontours=TRUE, pts=NULL,
          annot=annot, annot.cex=annot.cex, at=datarange, col.regions=cols,
          corners=p$corners, spatial.domain=p$spatial.domain )
        print(lp)
        dev.off()

        print( fn )
    }
    return( "Completed")
  }
}
