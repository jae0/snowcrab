

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields

# some issues running in MS Windows .. might need to run in Linux
#Virtual box install of Ubuntu or Debian is likely easiest option

year.assessment = 2020

# adjust based upon RAM requirements and ncores
require(INLA)
inla.setOption(num.threads= floor( parallel::detectCores() / 2) )
inla.setOption(blas.num.threads= 2 )


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)
 p = bio.snowcrab::snowcrab_parameters( project_class="carstm", assessment.years=1999:year.assessment )

  # p$modeldir = "..."  # use this to specifiy alt location to save model output files



# ------------------------------------------------
# Part 2 -- polygon structure
  if (0) {
    # create if not yet made
    for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
    xydata = snowcrab.db( p=p, DS="areal_units_input", redo=TRUE )
    sppoly = areal_units( p=p, redo=TRUE )  # create constrained polygons with neighbourhood as an attribute
    coastLayout = aegis.coastline::coastline_layout( p=p, redo=TRUE )
    MS = NULL
  }
  sppoly = areal_units( p=p )  # to reload
  # plot(sppoly)
  # spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout )
  dev.new(); 
  plot( sppoly[, "au_sa_km2"], main="AUID", breaks = "quantile", nbreaks = 8, pal=sf.colors(8) )

 
# -------------------------------------------------
# Part 8 -- Snow crab abundance -- main mode used for production
  M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  #To compare values of M, run the following line:
  #load(paste(p$modeldir, "M.summary.rdata", sep="/"))

  M = NULL; gc()


  fit = carstm_model( p=p, M='snowcrab.db( p=p, DS="carstm_inputs" )' ) # 151 configs and long optim .. 19 hrs
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


  require(aegis.coastline)
  coastline = coastline_db( p=p, DS="eastcoast_gadm" )
  coastline = st_transform( coastline, st_crs(p$aegis_proj4string_planar_km) )

  # depth contours
  require(aegis.polygons)
  isobaths = aegis.bathymetry::isobath_db( p=p, depths=c(50, 100, 200, 400, 800)  )
  isobaths = st_transform( isobaths, st_crs(p$aegis_proj4string_planar_km) )

  time_match = list( year=as.character(2020)  )

  vn = paste(p$variabletomodel, "predicted", sep=".")
     carstm_map(  res=res, vn=vn, time_match=time_match , 
          coastline=coastline,
          isobaths=isobaths,
          main=paste("Predicted abundance", paste0(time_match, collapse="-") )  
     )
    )


  # map all :
  vn = paste(p$variabletomodel, "predicted", sep=".")
  outputdir = file.path( gsub( ".rdata", "", dirname(res$fn_res) ), "figures", vn )
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )


  for (y in res$yrs ){

      time_match = list( year=as.character(y)  )
      fn_root = paste( "Bottom temperature",  paste0(time_match, collapse=" - ") )
      fn = file.path( outdir, paste(fn_root, "png", sep=".") )
      o = carstm_map(  res=res, vn=vn, time_match=time_match, 
#          breaks = seq(-0.5, 0.5, by=0.1),
          coastline=coastline,
          isobaths=isobaths,
          main=paste("Predicted abundance", paste0(time_match, collapse="-") )  
       main=fn_root 
      )
      png( filename=fn, width=3072, height=2304, pointsize=40, res=300  ); print(o); dev.off()
  }
  

  RES = snowcrab.db(p=p, DS="carstm_output_compute" )
  # RES = snowcrab.db(p=p, DS="carstm_output_timeseries" )

  bio = snowcrab.db(p=p, DS="carstm_output_spacetime_biomass" )
  num = snowcrab.db(p=p, DS="carstm_output_spacetime_number" )

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
  p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")
    sppoly = areal_units( p=p )  # to reload


  # map it ..mean density
  vn = "pred"
  slot(sppoly, "data")[,vn] = bio[,"2019"]
  brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
  spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )

  #to save plots of the last six years:
  #Needs to be run in Windows or outside Linux Rstudio for now
  recent=as.character((year.assessment-6): year.assessment)

  plot.dir=paste(p$modeldir, "prediction.plots", year.assessment, sep="/" )

  for (x in recent){
    slot(sppoly, "data")[,vn] = bio[,x]
    brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
    o=spplot( sppoly, vn, col.regions=p$mypalette, main=x, at=brks, sp.layout=p$coastLayout, col="transparent" )
    fn=paste(x,"biomass",  "pdf", sep=".")
    outfile=paste(plot.dir, fn, sep="/")
    pdf(outfile)
    print(o)
    dev.off()
  }


  plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
  plot( fit$marginals.hyperpar$"Precision for auid", type="l")
  plot( fit$marginals.hyperpar$"Precision for setno", type="l")






# end
