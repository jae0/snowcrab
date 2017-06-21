
  
  if (!exists("current.year")) current.year=year(Sys.Date())
  p = bio.snowcrab::load.environment( year.assessment=current.year)

# ------------------------------------------
# Time-series: all interpolated data estimated from interpolated analysis
  figure.interpolated.results( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ), alt.zero.y=T )

# ------------------------------------------
# Time-series: immature male numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.immature.male( p, outdir=file.path( p$annual.results, "timeseries", "interpolated" ) )

 # ------------------------------------------
# Time-series: immature male skip-moulter numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.immature.skipmoulter.male( p, outdir=file.path( p$annual.results,  "timeseries", "interpolated" ) )

# ------------------------------------------
# Time-series: mature CC12 male numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.mature.CC12.male( p, outdir=file.path( p$annual.results,  "timeseries", "interpolated" ) )

# ------------------------------------------
# Time-series: mature CC34 male numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.mature.CC34.male( p, outdir=file.path( p$annual.results, "timeseries", "interpolated" ) )

# ------------------------------------------
# Time-series: mature CC5 male numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.mature.CC5.male( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ) )


# ------------------------------------------
# Time-series: Exploitation rate
  figure.exploitation.rate( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ) )
  figure.exploitation.rate( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ), CFA4X.exclude.year.assessment=F )

# ------------------------------------------
# Time-series: Generic plots of all interpolated data
  # figure.timeseries.interpolated( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ) )

# ------------------------------------------
# Time-series: Habitat variations (surface area of snow crab habitat)

# p$model.type = "gam.full" # choose method for habitat model :
# p$model.type = "gam.simple" # choose method for habitat model :  no longer used

p$habitat.threshold.quantile = 0.05 # quantile at which to consider zero-valued abundance
p$prediction.dyear = 9/12 # predict for ~ Sept 1
figure.timeseries.snowcrab.habitat( p=p)

# ------------------------------------------
# Time-series: Habitat variations (mean temperature of snow crab habitat)
figure.timeseries.snowcrab.habitat.temperatures( p=p)

# ------------------------------------------
# Fecundity estimated indirectly via total number of females of primi and muli and applying mean egg production from allometric estimate
  fecundity.indirect( p=p,  outdir=file.path( p$annual.results, "timeseries", "survey" ) )

# ------------------------------------------
# Fecundity estimated directly via interpolation and individual-based fecundity estimate
  figure.timeseries.fecundity( p=p, outdir=file.path( p$annual.results, "timeseries",  "interpolated"  ) )




