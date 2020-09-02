

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields

# some issues running in MS Windows .. might need to run in Linux
#Virtual box install of Ubuntu or Debian is likely easiest option

year.assessment = 2019


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)
 p = bio.snowcrab::snowcrab_carstm( DS="parameters", assessment.years=1999:year.assessment )

  # misc run params adjustments here:
  p$inla_num.threads = 6
  p$inla_blas.num.threads = 6

  plot.dir=paste(p$modeldir,"prediction.plots", year.assessment, sep="/" )

# ------------------------------------------------
# Part 2 -- polygon structure
  if (0) {
    # create if not yet made
    for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygons_managementarea( species="snowcrab", au))
    sppoly = areal_units( p=p, redo=TRUE )  # create constrained polygons with neighbourhood as an attribute
    coastLayout = aegis.coastline::coastline_layout( p=p, redo=TRUE )
    MS = NULL
  }
  sppoly = areal_units( p=p )  # to reload
  # plot(sppoly)
  # spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout )
  x11(); spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout,  col.regions=RColorBrewer::brewer.pal(8, "Accent") )


# -------------------------------------------------
# Part 3 -- create covariate field for bathymetry
# bathymetry -- ensure the data assimilation in bathymetry is first completed :: 01.bathymetry_data.R
  pB = bathymetry_carstm( p=p, DS="parameters", variabletomodel="z" )
  M = bathymetry.db( p=pB, DS="aggregated_data" , redo=TRUE )
  M = bathymetry_carstm( p=pB, DS="carstm_inputs", redo=TRUE  ) # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=pB, M='bathymetry_carstm( p=pB, DS="carstm_inputs" )', DS="redo"  ) # run model and obtain


# -------------------------------------------------
# Part 4 -- create covariate field for  substrate
# ensure the data assimilation in substrate is first completed :: 01.substrate_data.R
  pS = substrate_carstm( p=p, DS="parameters", variabletomodel="substrate.grainsize" )
  M = substrate.db( p=pS, DS="aggregated_data", redo=TRUE )  # used for data matching/lookup in other aegis projects that use substrate
  M = substrate_carstm( p=pS, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=pS, M='substrate_carstm( p=pS, DS="carstm_inputs")', DS="redo" )  # run model and obtain predictions


# -------------------------------------------------
# Part 5 -- create covariate field for temperature
# ensure the data assimilation in temperature is first completed :: 01.temperature_data.R
  pT = temperature_carstm( p=p, DS="parameters", variabletomodel="t" )
  M = temperature.db( p=pT, DS="aggregated_data", redo=TRUE )  #  used for data matching/lookup in other aegis projects that use temperature
  M = temperature_carstm( p=pT, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=pT, M='temperature_carstm( p=pT, DS="carstm_inputs")', DS="redo"  ) # run model and obtain predictions

  # to plot predicted temperature maps for last six years

  res = carstm_summary( p=pT )
  recent=as.character((year.assessment-6): year.assessment)
  vn = paste(pT$variabletomodel, "predicted", sep=".")
  for (x in recent){
    fn=paste(x,"t",  "pdf", sep=".")
    outfile=paste(plot.dir, fn, sep="/")
    each.plot=   carstm_plot( p=pT, res=res, vn=vn, time_match=list(year=x, dyear="0.85" ) )
    pdf(outfile)
    print(each.plot)
    dev.off()
  }


# -------------------------------------------------
# Part 6 -- create covariate field for species composition 1
# ensure that survey data is assimilated : bio.snowcrab::01snowcb_data.R, aegis.survey::01.surveys.data.R , etc.
  pPC1 = speciescomposition_carstm( p=p, DS="parameters", variabletomodel="pca1" )
  M = speciescomposition_carstm( p=pPC1, DS="carstm_inputs",  redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=pPC1, M='speciescomposition_carstm( p=p, DS="carstm_inputs" )', DS="redo"   ) # run model and obtain predictions


# -------------------------------------------------
# Part 7 -- create covariate field for species composition 2
# ensure that survey data is assimilated : bio.snowcrab::01snowcb_data.R, aegis.survey::01.surveys.data.R ,
  pPC2 = speciescomposition_carstm( p=p, DS="parameters", variabletomodel="pca2")
  M = speciescomposition_carstm( p=pPC2, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=pPC2, M='speciescomposition_carstm( p=p, DS="carstm_inputs" )', DS="redo"  ) # run model and obtain predictions


# -------------------------------------------------
# Part 8 -- Snow crab abundance -- main mode used for production
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  #To compare values of M, run the following line:
  #load(paste(p$modeldir, "M.summary.rdata", sep="/"))

  M = NULL; gc()

  fit = carstm_model( p=p, M='snowcrab_carstm( p=p, DS="carstm_inputs" )' ) # 151 configs and long optim .. 19 hrs
  # fit = carstm_model( p=p, DS="carstm_modelled_fit")

  summary(fit)
  res = carstm_summary( p=p )

  vn = paste(p$variabletomodel, "predicted", sep=".")
  carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000" ) )     # maps of some of the results

  plot(fit)
  plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )
  s = summary(fit)
  s$dic$dic
  s$dic$p.eff

  RES = snowcrab_carstm(p=p, DS="carstm_output_compute" )
  # RES = snowcrab_carstm(p=p, DS="carstm_output_timeseries" )

  bio = snowcrab_carstm(p=p, DS="carstm_output_spacetime_biomass" )
  num = snowcrab_carstm(p=p, DS="carstm_output_spacetime_number" )

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
  sppoly@data[,vn] = bio[,"2019"]
  brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
  spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )

  #to save plots of the last six years:
  #Needs to be run in Windows or outside Linux Rstudio for now
  recent=as.character((year.assessment-6): year.assessment)

  for (x in recent){
    sppoly@data[,vn] = bio[,x]
    brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
    each.plot=spplot( sppoly, vn, col.regions=p$mypalette, main=x, at=brks, sp.layout=p$coastLayout, col="transparent" )
    fn=paste(x,"biomass",  "pdf", sep=".")
    outfile=paste(plot.dir, fn, sep="/")
    pdf(outfile)
    print(each.plot)
    dev.off()
  }


  plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
  plot( fit$marginals.hyperpar$"Precision for auid", type="l")
  plot( fit$marginals.hyperpar$"Precision for setno", type="l")






# end
