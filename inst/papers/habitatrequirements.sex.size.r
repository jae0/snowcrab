
  # basic processing to compare sex and size-dependent habitat requirements of snow crab

  # define categories of interest


  p = bio.snowcrab::load.environment.r()

  p$regions ="cfaall"
  p$vars.to.model = c("male.large.mass", "male.small.mass", "female.large.mass", "female.small.mass" )

  p$do.parallel =F
  p$clusters = rep("localhost",8)
  # p$clusters = rep("tethys", 7 )
  # p$clusters = rep("kaos",23 )


      # Define controlling parameters

      p$model.type = "gam.full" # choose method for habitat model :
      p$habitat.threshold.quantile = 0.05 # quantile at which to consider zero-valued abundance
      p$optimizers = c( "nlm", "perf", "bam", "bfgs", "newton", "Nelder-Mead" )  # used by GAM
			p$prediction.dyear = lubridate::decimal_date( lubridate::ymd("0000/Sep/01")) # predict for ~ Sept 1
      p$threshold.distance = 14  # limit to extrapolation/interpolation in km


      # ---------------------
      # model habitat and intermediate predictions
      # ~ 1 MB / process  .. use all 24 CPUs
      # Parameterize habitat space models for various classes,
      # see: habitat.model.db( DS="habitat.redo" ... )

      # p$clusters = rep( "localhost", 1)
      # p$clusters = c( rep( "nyx.beowulf", 24), rep("tartarus.beowulf", 24), rep("kaos", 24 ) )
      p = make.list( list(v=p$vars.to.model, yrs=p$yrs ), Y=p )
      parallel.run( habitat.model.db, DS="habitat.redo", p=p )

      # or
      # habitat.model.db( DS="habitat.redo", p=p )





      # ---------------------
      # habitat surface area estimation for R0.mass from 1970 to present --- for timeseries and direct prediction maps
      p$clusters = c( rep( "localhost", 20) )
      # p$clusters = rep( "localhost", 2)

      p = make.list( list(v=c("R0.mass.environmentals.only", "R0.mass") ), Y=p )
        habitat.model.db (DS="habitat.redo", p=p )


