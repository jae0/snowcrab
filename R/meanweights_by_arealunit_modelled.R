

meanweights_by_arealunit_modelled = function( p=NULL, redo=FALSE, returntype="predictions_only" ) {
  # find mean weight for each stratum and year .. fill where possible if required
 
  if (redo) {
    message( "Re-creating mean size estimates from data .. this will take about 2 hrs ...")
    # must be fetched beofre modfying the model_label and variable name
    M = snowcrab.db( p=p, DS="carstm_inputs" ) 
    M$meansize  = M$totwgt / M$totno  # note, these are constrained by filters in size, sex, mat, etc. .. in the initial call
  }

  # alter model formula / storage location
  p_mw = p
  
  p_mw$variabletomodel = "meansize"
  p_mw$carstm_model_label = paste( p_mw$variabletomodel, p_mw$areal_units_type, p_mw$selection$type, sep="_")  

  p_mw$carstm_model_formula = as.formula( paste(
      p_mw$variabletomodel, ' ~ 1',
            ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
            ' + f( yr, model="ar1",  hyper=H$ar1 ) ',
            ' + f( inla.group( t, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( inla.group( z, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( inla.group( substrate.grainsize, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( inla.group( pca1, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( inla.group( pca2, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2 ) ',
            ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), group=time_space, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ) )

  p_mw$carstm_model_family =  "gaussian"   


  if (redo) {
 
    fit = carstm_model( p=p_mw, M=M ) # 151 configs and long optim .. 19 hrs
  } else{
    fit = carstm_model( p=p_mw, DS="carstm_modelled_fit" )  # extract currently saved model fit
  }
    
  res = carstm_model( p=p_mw, DS="carstm_modelled_summary"  ) # to load currently saved results


    if (0) {

          plot_crs = p_mw$aegis_proj4string_planar_km
          coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs )
          isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs )
          managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
        
          time_match = list( year=as.character(2020)  )

          carstm_map(  res=res, 
              vn=paste(p_mw$variabletomodel, "predicted", sep="."), 
              time_match=time_match, 
              coastline=coastline,
              managementlines=managementlines,
              isobaths=isobaths,
              main=paste("Predicted numerical abundance", paste0(time_match, collapse="-") )  
          )
            
    }
  if (returntype == "all" ) return( res )

  return(res$meansize.predicted)
}

