

# snow crab using alt AUs
# using global work flow : each data layer has its own optimal solution/au
# it is therefore a hybrid of sorts, using stmv in high data density situations (bathymetry)
# and data-dependent au's for other data layers

# -- results can be reused outside of workflow


# -------------------------------------------------
# 1. set up polygon parameters

  yrs = 1999:2019

  #  areal_units_type = "lattice"
  areal_units_type = "tesselation"

  # 5 - 10 works well .. mean and variances stabilize
  # 3 is too low for ts analysis,...  esp for temperature
  # areal_units_constraint_nmin  = floor(length(yrs) / 3) # = 6
  # areal_units_constraint_nmin = 3
  # areal_units_constraint_nmin = length(yrs)
  # areal_units_constraint_nmin = 4
  # areal_units_constraint_nmin = 10
  # areal_units_constraint_nmin = 15
  # areal_units_constraint_nmin = 30
  areal_units_constraint_nmin = 10

  #   areal_units_resolution_km = 2
  #   areal_units_resolution_km = 5
  #   areal_units_resolution_km = 10 # too disconnected (at nmin=3)
  #   areal_units_resolution_km = 15 # well connected at nmin=3,
  #  areal_units_resolution_km =  20
  #   areal_units_resolution_km =  25  # default snow crab
  #   areal_units_resolution_km =  30
  #   areal_units_resolution_km =  40
  areal_units_resolution_km = 1

  p = bio.snowcrab::snowcrab_parameters(
    DS="carstm",
    assessment.years=1999:2019,
    modeldir = project.datadirectory("bio.snowcrab", "modelled", "testing" ),  ## <--- important: specify save location
    carstm_model_label = paste( "testing", areal_units_type, areal_units_resolution_km, areal_units_constraint_nmin, sep="_" ),
    inputdata_spatial_discretization_planar_km = 1,
    boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)), # bounding box for plots using spplot
    areal_units_proj4string_planar_km = projection_proj4string("utm20"), # set up default map projection
    areal_units_constraint = "snowcrab",
    areal_units_constraint_nmin = areal_units_constraint_nmin,
    areal_units_type= areal_units_type,
    areal_units_resolution_km = areal_units_resolution_km,
    sa_threshold_km2 = 5
  )


  (p$carstm_model_label ) # required to use for abundance estimation
  # [1] carstm_model_label="testing_snowcrab_polygons_tesselation_1_3"

  if (0) coastLayout = aegis.coastline::coastline_layout( p=p, redo=TRUE )

  p = c(p, aegis.coastline::coastline_layout( p=p ) )
  # p$mypalette = RColorBrewer::brewer.pal(9, "YlOrRd")
  p$mypalette = RColorBrewer::brewer.pal(9, "Spectral")[9:1]


  # p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")   # oblique mercator, centred on Scotian Shelf rotated by 325 degrees



  # create polygons
  if (0) {
    # ensure if polys exist and create if required
    for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) {
      plot(polygon_managementareas( species="snowcrab", au))
    }
  }

  sppoly = areal_units( p=p )


  if (0) {
    # verify that params are correct
    sppoly = areal_units( p=p, redo=TRUE )  # to create

    plot(sppoly[,"AUID"], col="orange")
    plot( slot(sppoly, "nb"), coords=st_centroid(st_geometry( as(sppoly, "sf")) ),  col="green", add=T )

    dev.new(); spplot( sppoly, "au_sa_km2", main="SA_km2", sp.layout=p$coastLayout )
    dev.new(); spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout,  col.regions=RColorBrewer::brewer.pal(8, "Accent") )
    length(sppoly)
    summary(sppoly$au_sa_km2)

  }


# -------------------------------------------------
# 2. run global interpolations in each layer

# -- make sure each of the following have been run, where necesssary
# bathymetry.R , etc



# -------------------------------------------------
# 3  collect aegis data by areal units from global sourcesSnow and add crab anbundance --


if (0) {
  # use alternate model -- zero-inflated1

    p$carstm_model_label = "zeroinflated"  #unique to this project ... to permit alt model forms/variations within the same overall carstm

     p$carstm_model_formula = as.formula( paste(
      p$variabletomodel, ' ~ 1',
        ' + offset( log(data_offset))',
        ' + f( dyri, model="ar1", hyper=H$ar1 )',
        ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
        ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
        ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
        ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
        ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
        ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group))'
    ) )
    p$carstm_model_family = "zeroinflatedpoisson1"

}






  M = snowcrab_db( p=p, DS="carstm_inputs_hybrid", redo=TRUE )  # will redo if not found
  fit = carstm_model( p=p, M='snowcrab_db( p=p, DS="carstm_inputs_hybrid" )' ) # 151 configs and long optim .. 19 hrs
  RES = snowcrab_db(p=p, DS="carstm_output_compute" )

  if (0) {

      fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
      summary(fit)

      res = carstm_model( p=p, DS="carstm_modelled_summary" )

      vn = paste(p$variabletomodel, "predicted", sep=".")
      carstm_map( res=res, vn=vn, time_match=list(year="2019" ), 
        # at=seq(-2, 10, by=2),          
        sp.layout = p$coastLayout, 
        col.regions = p$mypalette, 
        main=paste( vn, paste0(time_match, collapse="-") )  
      )
     
     
      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )
      s = summary(fit)
      s$dic$dic
      s$dic$p.eff

      # maps of some of the results
      time_match = NULL
      time_match = list(year="2000")
      time_match = list(year="2000", dyear="0.65" )
      
      vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
      carstm_map( res=res, vn=vn, time_match=time_match, 
        # at=seq(-2, 10, by=2),          
        sp.layout = p$coastLayout, 
        col.regions = p$mypalette, 
        main=paste( vn, paste0(time_match, collapse="-") )  
      )


      vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
      carstm_map( res=res, vn=vn, time_match=time_match, 
        # at=seq(-2, 10, by=2),          
        sp.layout = p$coastLayout, 
        col.regions = p$mypalette, 
        main=paste( vn, paste0(time_match, collapse="-") )  
      )


      vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
      carstm_map( res=res, vn=vn, time_match=time_match, 
        # at=seq(-2, 10, by=2),          
        sp.layout = p$coastLayout, 
        col.regions = p$mypalette, 
        main=paste( vn, paste0(time_match, collapse="-") )  
      )
     
      
      RES = snowcrab_db(p=p, DS="carstm_output_timeseries" )
      bio = snowcrab_db(p=p, DS="carstm_output_spacetime_biomass" )
      num = snowcrab_db(p=p, DS="carstm_output_spacetime_number" )


      plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
      lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

      plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
      lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

      plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
      lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

      plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
      lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )


      p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot

      p$coastLayout = aegis.coastline::coastline_layout(p=p)
      # p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")

      sppoly = areal_units( p=p )  # to reload

      # map it ..mean density
      vn = "pred"
      slot(sppoly, "data")[,vn] = bio[,"2019"]
      brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
      spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


      plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
      plot( fit$marginals.hyperpar$"Precision for auid", type="l")
      plot( fit$marginals.hyperpar$"Precision for setno", type="l")
  }


}
}



# ----

# estimate numbers for "recruitment"

  p$selection$type = "number"

  p$selection$biologicals = list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    sex=0, # male
    # mat=1, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 80, 94.99 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )


  M = snowcrab_db( p=p, DS="carstm_inputs_hybrid", redo=TRUE )  # will redo if not found
 
  fit = carstm_model( p=p, M='snowcrab_db( p=p, DS="carstm_inputs_hybrid" )' ) # 151 configs and long optim .. 19 hrs
  # fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  
  
  res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
  res$summary$dic$dic
  res$summary$dic$p.eff
  res$dyear

RES = snowcrab_db(p=p, DS="carstm_output_compute" )



# for the abundance run :


----


# runtype="biomass"
runtype="number"
# runtype="presence_absence"

# ---------------------

# basic spatial parameters
p = carstm::carstm_parameters(
  id ="Snow crab",
  speciesname = "Snow crab",
  groundfish_species_code = 2526,
  spatial_domain="snowcrab",
  yrs = 1999:2018,
  areal_units_resolution_km = 20, # km
  areal_units_type = "lattice",
  # areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia"),   # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
  areal_units_proj4string_planar_km = projection_proj4string("utm20"),
  boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)), # bounding box for plots using spplot
  trawlable_units = "towdistance",  # <<<<<<<<<<<<<<<<<< also:  "standardtow", "sweptarea" (for groundfish surveys)
  libs = RLibrary ( "sp", "spdep", "rgeos", "INLA", "raster", "aegis", "aegis.polygons", "aegis.coastline", "aegis.survey", "bio.taxonomy", "carstm" )
)

p = aegis.survey::survey_parameters(
  p=p,
  selection=list(
    type = runtype,
    biologicals=list(
      spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
      sex=0, # male
      mat=1, # do not use maturity status in groundfish data as it is suspect ..
      len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
      ranged_data="len"
    ),
    survey=list(
      data.source = ifelse (runtype=="number", c("snowcrab"), c("snowcrab", "groundfish")),
      yr = p$yrs,      # time frame for comparison specified above
      settype = 1, # same as geartype in groundfish_survey_db
      polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
      strata_toremove = NULL,  # emphasize that all data enters analysis initially ..
      ranged_data = c("dyear")  # not used .. just to show how to use range_data
    )
  )
)

p$lookupvars = c("t", "tsd", "tmax", "tmin", "degreedays", "z",  "dZ", "ddZ", "substrate.grainsize" )



# --------------------------------
# Get the data

  set = bio.snowcrab::snowcrab_stmv( p=p, DS="input_data", alldata=TRUE, redo=reset_input_data )
  set$Y = set$totno  # unadjusted value is used as we are usinmg offsets ...
  set$data_offset  = 1 / set[, ifelse( p$selection$type=="number", "cf_set_no", "cf_set_wgt")]  # as "sa"
  set$data_offset[which(!is.finite(set$data_offset))] = median(set$data_offset, na.rm=TRUE )  # just in case missing data
  set$tag = "observations"


  # extract covariates and supplent survey data via lookups
  # in snow crab this is already done .. but for other edatastreams this would be useful here
  if (0) {
    # if needed
    set = aegis_db_lookup(
      X=set,
      lookupvars=p$lookupvars,
      xy_vars=c("lon", "lat"),
      time_var="timestamp"
    )
  }





## ----------------------------------
# covariate estimates for prediction in strata and year
# collapse PS vars with time into APS (and regrid via raster)
  APS = aegis_db_extract(
    vars=p$lookupvars,
    yrs=p$yrs,
    spatial_domain=p$spatial_domain,
    dyear=p$prediction_dyear,
    areal_units_resolution_km=p$areal_units_resolution_km,
    aegis_proj4string_planar_km=p$aegis_proj4string_planar_km,
    returntype="data.frame",
    redo = reset_input_data
  )

  APS$yr = as.numeric( APS$year)
  APS$Y = NA
  APS$data_offset = 1  # force to be density n/km^2
  APS$tag = "predictions"


  # AUID reset to be consistent in both data and prediction areal units
  crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
  sppoly = st_transform(sppoly, crs=crs_lonlat )

  set$AUID = st_points_in_polygons(
    pts = st_as_sf( set, coords=c("lon","lat"), crs=crs_lonlat ),
    polys = sppoly[, "AUID"],
    varname="AUID"
  )


  APS = planar2lonlat(APS, p$aegis_proj4string_planar_km)

  APS$AUID = st_points_in_polygons(
    pts = st_as_sf( APS, coords=c("lon","lat"), crs=crs_lonlat ),
    polys = sppoly[, "AUID"],
    varname="AUID"
  )

  #  good data
  ok = which(
    is.finite(set[,p$variabletomodel]) &   # INLA can impute Y-data
    is.finite(set$data_offset) &
    !is.na(set$AUID)
  )


# construct meanweights matrix
weight_year = meanweights_by_arealunit( set=set, AUID=as.character( sppoly$AUID ), yrs=p$yrs, fillall=TRUE, annual_breakdown=TRUE )
# weight_year = weight_year[, match(as.character(p$yrs), colnames(weight_year) )]
# weight_year = weight_year[ match(as.character(sppoly$AUID), rownames(weight_year) )]


varstokeep = unique( c( "Y", "AUID", "yr", "data_offset", "tag", p$lookupvars) )

M = rbind( set[ok, varstokeep], APS[,varstokeep] )

M = M[ which(
      is.finite(M$data_offset) &
      !is.na(M$AUID)
    ) , ]

M$yr_factor = factor( as.character(M$yr) )
# M$AUID  = factor( as.character(M$AUID), levels=levels( sppoly$AUID ) )

region.id = slot( slot(sppoly, "nb"), "region.id" )
M$auid = match( M$AUID, region.id )

M$AUID_character = M$AUID
M$AUID = set$auid  -- temp fix to get numbers to match bym geometry

M$strata  = as.numeric( M$AUID)
M$year  = as.numeric( M$yr_factor)
M$iid_error = 1:nrow(M) # for inla indexing for set level variation


M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$tsd[!is.finite(M$tsd)] = median(M$tsd, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$tmin[!is.finite(M$tmin)] = median(M$tmin, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$tmax[!is.finite(M$tmax)] = median(M$tmax, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$degreedays[!is.finite(M$degreedays)] = median(M$degreedays, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$dZ[!is.finite(M$dZ)] = median(M$dZ, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$ddZ[!is.finite(M$ddZ)] = median(M$ddZ, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$substrate.grainsize[!is.finite(M$substrate.grainsize)] = median(M$substrate.grainsize, na.rm=TRUE )  # missing data .. quick fix .. do something better


M$ti = discretize_data( M$t, p$discretization$t )
M$tisd = discretize_data( M$tsd, p$discretization$tsd )
M$timin = discretize_data( M$tmin, p$discretization$tmin )
M$timax = discretize_data( M$tmax, p$discretization$tmax )
M$di = discretize_data( M$t, p$discretization$degreedays )
M$zi = discretize_data( M$t, p$discretization$z )
M$zid = discretize_data( M$t, p$discretization$dZ )
M$zidd = discretize_data( M$t, p$discretization$ddZ )
M$si = discretize_data( M$t, p$discretization$substrate.grainsize )


totest = setdiff(1:ncol(M), which(names(M) %in% c("Y", "AUID", "tag", "yr_factor") ))
ii = which(is.finite(rowSums(M[,totest])))

M = M[ii,]

# ---------------------
# generic PC priors
m = log( {set$Y / set$data_offset}[ok] )
m[!is.finite(m)] = min(m[is.finite(m)])






