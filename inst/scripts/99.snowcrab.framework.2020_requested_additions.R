
## Plot maps of residuals of numbers per set obs vs pred


year.assessment = 2018
#To add a title to any carstm_plot, please see below example
#carstm_plot( p=p, res=res, vn=vn, main=list(label="my plot title", cex=2) )


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)

 p = bio.snowcrab::snowcrab_carstm( DS="parameters", assessment.years=1999:year.assessment )

  # misc run params adjustments here:
  p$inla_num.threads = 6
  p$inla_blas.num.threads = 6

  plot.dir=paste(p$modeldir,"prediction.plots", year.assessment, sep="/" )



# extract results and examine
  fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  summary(fit)

  res = carstm_summary( p=p, operation="compute"  ) #
  res = carstm_summary( p=p, operation="load"  )



    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found
    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))


    # do this immediately to reduce storage for sppoly (before adding other variables)
    M = snowcrab.db( p=p, DS="biological_data" )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry

    #January 2020 samples create a problem- 2019 survey, 2020 year
    #Shift these samples back to late december by removing 16 days
    i=which((lubridate::year (M$timestamp)==2020) & (lubridate::month(M$timestamp)==1))
    #M$tiyr[i]=M$tiyr[i]-(16/365.25)
    M$timestamp[i]=M$timestamp[i]-1382400
    M$tiyr=lubridate::decimal_date(M$timestamp)

    # M$totno = M$totno_adjusted / M$cf_set_no   # convert density to counts
    # M$totwgt = M$totwgt_adjusted / M$cf_set_mass # convert density to total wgt

    # M$data_offset = 1 / M$cf_set_no  ## offset only used in poisson model


    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")

    M$AUID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
    M = M[!is.na(M$AUID),]

    names(M)[which(names(M)=="yr") ] = "year"
    # M = M[ which(M$year %in% p$yrs), ]
    # M$tiyr = lubridate::decimal_date ( M$timestamp )
    # M$dyear = M$tiyr - M$year

    MM = snowcrab_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found


    obsMM = MM[MM$tag=="observations",]
    plocs = MM[MM$tag=="predictions",]

    obs = M
    obs$density = obs$totno / obs$data_offset

    rr = as.data.frame.table(res$totno.predicted)
    rr$AUID = as.character( rr$AUID)
    rr$year = as.numeric( as.character( rr$year) )

    obs = merge( obs, rr, by=c("year", "AUID"), all.x=TRUE, all.y=FALSE )
    obs$no_resid = obs$density - obs$Freq
    obs$no_resid_per_set = obs$no_resid * obs$data_offset
    obs$yr = obs$year


    vn = "no_resid"
    vn = "no_resid_per_set"

    resol = p$pres

    B = bathymetry.db(p=p, DS="baseline")  # 1 km (p$pres )

    for ( y in  2000:2018 ) {
        ii = which( obs$yr==y)
        if ( length(ii) > 3 ) {
        fn = file.path( p$project.outputdir, "residuals", y, "rdata", sep=".")
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = map_simple( toplot=obs[ ii, c("plon","plat", vn) ], plotarea=B, resol=1, theta=7.5, filterdistances=7.5, vn=vn, annot=paste("Residuals", y) )
        print(lp)
        dev.off()
      }
    }





