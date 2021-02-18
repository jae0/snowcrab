

# Snow crab --- Areal unit modelling of habitat   

# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study


  year.assessment = 2020

  p = bio.snowcrab::snowcrab_parameters( 
    project_class="carstm", 
    assessment.years=2000:year.assessment,  
    areal_units_type="tesselation"
  )

  p$selection$type = "presence_absence"
  p$selection$biologicals=list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    sex=0, # male
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )
  p$selection$survey=list(
    data.source = c("snowcrab", "groundfish", "logbook"),
    yr = p$assessment.years,      # time frame for comparison specified above
    settype = 1, # same as geartype in groundfish_survey_db
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
    # ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
  p$variabletomodel = "pa"
  
  p$carstm_model_label = "nonseparable_space-time_pa_fishable_zeroinflatedbinomial1"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ) )

  p$carstm_model_family = "binomial"  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
  p$carstm_model_inla_control_familiy = list(control.link=list(model='logit'))

#  p$carstm_model_family  = "zeroinflatedbinomial1", #  "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
#  p$carstm_model_inla_control_familiy = NULL



  sppoly = areal_units( p=p )  # to reload
  plot( sppoly[, "au_sa_km2"]  )


  M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=p, M=p$modeldata )




  # fit = carstm_model( p=p, DS="carstm_modelled_fit")

    # extract results
    if (0) {
      # very large files .. slow 
      fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    }


  res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
  res$summary$dic$dic
  res$summary$dic$p.eff
  res$dyear


  plot_crs = p$aegis_proj4string_planar_km
  coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs )
  isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs )
  managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
 
  carstm_map(  res=res, 
      vn=paste(p$variabletomodel, "predicted", sep="."), 
      time_match=list( year=as.character(2020)  ) , 
      coastline=coastline,
      managementlines=managementlines,
      isobaths=isobaths,
      main=paste("Habitat probability - mature male ", paste0(time_match, collapse="-") )  
  )
    

  # map all :
  vn = paste(p$variabletomodel, "predicted", sep=".")
  outputdir = file.path( gsub( ".rdata", "", dirname(res$fn_res) ), "figures", vn )
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  brks = pretty(  res[[vn]]  )

  for (y in res$year ){

      time_match = list( year=as.character(y)  )
      fn_root = paste("Predicted_abundance", paste0(time_match, collapse="-"), sep="_")
      fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

        carstm_map(  res=res, vn=vn, time_match=time_match, 
          breaks =brks,
          palette="viridis",
          coastline=coastline,
          isobaths=isobaths,
          managementlines=managementlines,
          main=paste("Habitat probability - mature male ", paste0(time_match, collapse="-") ),
          outfilename=fn
        )  

  }
  


  snowcrab.db(p=p, DS="carstm_output_compute" )
  
  RES = snowcrab.db(p=p, DS="carstm_output_timeseries" )

  pa = snowcrab.db(p=p, DS="carstm_output_spacetime_pa"  )

# plots with 95% PI

    plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="", ylim=c(0.2,0.4))
    lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

    plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="", ylim=c(0,1))
    lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

    plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="", ylim=c(0,1))
    lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

    plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="", ylim=c(0,1))
    lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

 

  # map it ..mean density

  sppoly = areal_units( p=p )  # to reload

  plot_crs = p$aegis_proj4string_planar_km
  coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs )
  isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs  )
  managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
 
  vn = paste("pa", "predicted", sep=".")

  outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.biomass.densitites" )

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  brks = pretty(  pa[]   )

  for (i in 1:length(p$yrs) ){
    y = as.character( p$yrs[i] )
    sppoly[,vn] = pa[,y] 
    fn = file.path( outputdir , paste( "biomass", y, "png", sep=".") )

      carstm_map(  sppoly=sppoly, vn=vn,    
        breaks=brks, 
        coastline=coastline,
        isobaths=isobaths,
        managementlines=managementlines,
        palette="-viridis",
        main=paste("Predicted habitat probability", y ),  
        outfilename=fn
      )
  }

  plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
  plot( fit$marginals.hyperpar$"Precision for auid", type="l")
  plot( fit$marginals.hyperpar$"Precision for setno", type="l")



 
# end
