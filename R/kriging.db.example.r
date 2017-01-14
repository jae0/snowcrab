
  # this is not operating .. it is an exmple only of the rough logic of kriging

  kriging.db.example = function( ip=NULL, DS, p=NULL, yrs, overwrite.threshold=NULL) {


    if (DS=="example") {
      # just a dummy example

  # --------------------------
  # 0. Initialise work space

  RLibrary( "mgcv", "chron", "lattice"  )
	p = snowcrab::load.environment()

  # --------------------------
  # 1. Define some additional starting parameters for debugging
  #    choose various over-rides: these are initially defined in parameters.r

  debug = F
  if (debug) {
    # identify areas to interpolate
      p$regions = c("cfanorth")
      p$regions = c("cfa4x", "cfanorth","cfasouth" )
      p$regions = c("cfa4x", "cfanorth","cfasouth", "cfaall" )

    # identify variables to interpolate and create abundance estimates

      p$vars.to.model = c( "R0.no", "R1.no")
      p$vars.to.model = c("R0.mass")
      p$vars.to.model = c("R0.mass", "male.large.mass", "male.small.mass", "female.large.mass", "female.small.mass" )
      p$vars.to.model = c("R0.mass", "R0.no", "R1.mass", "R1.no", "R2.no", "R3.no", "R4.no", "R5p.no", "totmass.male.com.CC1",
        "totno.female.berried","totno.female.imm", "totno.female.mat", "totno.female.primiparous","totno.female.multiparous",
        "fecundity", "totno.male.com", "totno.male.mat", "totno.male.imm","dwarf.no", "totno.male.skip.moulter",
        "totno.male.com.CC1to2",  "totno.male.com.CC5",  "male.large.mass", "male.small.mass", "female.large.mass", "female.small.mass", "male.large.np", "male.small.no", "female.large.no", "female.small.no"  )

    # modify cluster requirements
      p$do.parallel =F
      p$clusters = rep("localhost",8)
      # p$clusters = rep("tethys", 7 )
      # p$clusters = rep("kaos",23 )
      # p$clusters = rep("nyx",24 )
      # p$clusters = rep("tartarus",24 )
  }



  #-----------------------------------------------------
  # 3. KRIGING estimation of abundance and mapping

    # First, predict habitat surface locations: takes about 2hr for 15 yrs ..
    # snowcrab.habitat.db.old ( DS="habitat.redo", model.type=p$model.type, vclass=vclass, pyears=p$yearswithTdata, p=p )
    # or

      p$habitat.threshold.quantile = 0.05  # quantile at which to consider zero-valued abundance
      p$habitat.templates =  "R0.mass"
      p$model.type = "gam.full" # choose method for habitat model :

      model.habitat( model.type="gam.full.redo", p=p, plotdata=F ) # only a single habitat type is used with this method

      p = make.list( list( yrs=p$yearswithTdata ), Y=p )
      parallel.run( snowcrab.habitat.db.old, DS="habitat.redo", model.type=p$model.type, vclass=vclass, pyears=p$yearswithTdata, p=p )


      # V prepare prediction surface data sets for kriging analysis and GAM
      # snowcrab.habitat.db.old( DS="PS.redo", p=p, yrs=p$years.to.model,  model.type=p$model.type, vclass=vclass)
      # or

      p = make.list( list( yrs=p$years.to.model ), Y=p )
      parallel.run( snowcrab.habitat.db.old, DS="PS.redo", p=p, yrs=p$years.to.model, model.type=p$model.type, vclass=vclass )


			# determine potential habitat stats in historical data ( 1950 to present ) for timeseries plots
      # --- is this necessary ??
      prediction.surface( p, DS="annual.redo", vclass=vclass )


      # debug
      # x = snowcrab.habitat.db.old(DS="habitat", model.type=p$model.type, vclass=vclass, pyears=yr )
      # x11(); levelplot( x$fit ~ plon+plat, indicators.db( DS="static", p=p), aspect="iso" )

      if ( point.kriging) {
        # this is to point krige and then map results

        p$kriging.type =  "point.kriging"
        p$use.local.variograms = F
        p$plot.variogram = F  # set to true only if master thread is on the same computer as the parallel threads

        p$knmax = 40
        p$vgm.dist = p$vgm.dist [ which( p$vgm.dist <= 100 ) ]
        p$kriging.distance.buffer = 10  # use additional bounding data to this extent to predict
        p$kriging.extrapolation.limit = 10  # limit extrap to within x km of a surevey station


        # define compact list of variable year combinations for parallel processing
        p = make.list( list(y=p$years.to.model, r="cfaall", v=p$vars.to.model ), Y=p )

        # kriging.db( DS="UK.redo", p=p)
        # or
        parallel.run( kriging.db, DS="UK.redo", p=p)

        # the following is completed within the kriging.db but can be re-run this way:
        # map.krige.lattice( M=p)
        # parallel.run( map.krige.lattice, M=p)

        # test plot
        # PS = kriging.db( DS="UK.point.PS", p=list(v="R0.mass", y=p$year.assessment, r="cfanorth")  )
        # levelplot( R0.mass.pred ~ plon+plat, PS, aspect="iso")
        #
        # plot with GMT
        # krige.map ( p)  # ~ 1 day?
        # gmt.cleanup() # clean up any stragglers

      }


      if ( abundance.estimation.via.kriging ) {

        p$kriging.type = "block.conditional.sims"
        p$use.local.variograms = T
        p$n.conditional.sims = 150

        p$knmax = 40
        p$vgm.dist = p$vgm.dist [ which( p$vgm.dist <= 100 ) ]
        p$kriging.distance.buffer = 10  # use additional bounding data to this extent to predict
        p$kriging.extrapolation.limit = 10  # limit extrap to within x km of a surevey station

        debug.4X = F
        if (debug.4X) {
          problem.years = c(2006, 2009 )
          p$years.to.model = setdiff( c(2004:p$year.assessment), problem.years)   # highly problematic years
          p$years.to.model = c(2006, 2008, 2009)  # highly problematic years
          p$regions = "cfa4x"
          p$n.conditional.sims = 150
          p$vgm.dist =  p$vgm.dist [ which( p$vgm.dist <= 80 ) ]
          p$use.local.variograms = F
          p$knmax = 40
          p$kriging.distance.buffer = 10  # use additional bounding data to this extent to predict
          p$kriging.extrapolation.limit = 10  # limit extrap to within x km of a survey station
          p$vars.to.model = "R0.mass"
        }


        # define compact list of variable year combinations for parallel processing
        p = make.list( list(y=p$years.to.model, r=p$regions, v=p$vars.to.model ), Y=p )

        # kriging.db( DS="UK.redo", p=p )
        # kriging.db( DS="UK.redo", p=p, overwrite.threshold=NULL ) # overwrite all data files
        # kriging.db( DS="UK.redo", p=p, overwrite.threshold=-1 ) # complete only if no data files exist
        # kriging.db( DS="UK.redo", p=p, overwrite.threshold=100 ) # overwrite if over a certain number of days
        parallel.run( kriging.db, DS="UK.redo", p=p, overwrite.threshold=NULL  )
        # parallel.run( kriging.db, DS="UK.redo", p=p, overwrite.threshold=200  )  # ~ 0.4 GB / process (GSTAT)
        # ~ 0.4 GB / process (GSTAT);;
        # overwrite.threshold is # days required to trigger an over-write
        # overwrite.threshold=NULL  is overwrite all


        # collect all results into a single file and return:
        K = kriging.db( DS="UK.conditional.simulation.K", p=p  )
        table.view( K )

        if ( exists( "compare.density.with.biomass" ) ) {

          s = snowcrab.db( DS ="set.biologicals" )
          s1 = as.data.frame.table ( tapply( s$R0.mass, list(cfa=s$cfa, yr=s$yr ), mean, na.rm=T ) )
          names(s1) = c("cfa", "yr", "R0.mass.mean.density" )
          K = kriging.db( DS="UK.conditional.simulation.K", p=p  )
          o = K[ which(K$vars=="R0.mass"),]
          o1 = as.data.frame.table ( tapply( o$total, list(cfa=o$region, yr=o$yr ), mean, na.rm=T ) )
          names(o1) = c("cfa", "yr", "R0.mass.total.kriged" )
          o2 = as.data.frame.table ( tapply( o$surfacearea, list(cfa=o$region, yr=o$yr ), mean, na.rm=T ) )
          names(o2) = c("cfa", "yr", "sa" )
          o3 = merge( o1, o2, by=c("cfa", "yr") )
          o4 = merge( o3, s1, by=c("cfa", "yr") )
          sa = as.data.frame.table( tapply( o4$sa, o4$cfa, mean, na.rm=T) )
          names( sa ) = c("cfa", "sa.mean" )
          o = merge( o4, sa, by="cfa" )
          o$R0.mass.total.fixed = o$R0.mass.mean.density * o$sa.mean
          o$R0.mass.total.variable = o$R0.mass.mean.density * o$sa
          plot( o$R0.mass.total.variable, o$R0.mass.total.fixed, pch=19, xlab="Abundance (variable area)", ylab="Abundance (fixed area)" )
         #  text(  o$R0.mass.total.variable, o$R0.mass.total.fixed, o$yr, pos=4)
          abline(0,1)
          Pr("png", dname="~/", fname="tmp")
          hist(o$R0.mass.total.fixed/o$R0.mass.total.variable)

        }

        xyplot( total /1000~ yr | region, K, subset=which(K$vars=="R0.mass"), type="b", scales=list(y="free") )
        xyplot( total /1000~ yr | region, K, subset=which(K$vars=="R1.no"), type="b", scales=list(y="free") )
        xyplot( total /1000~ yr | region, K, subset=which(K$vars=="fecundity"), type="b", scales=list(y="free") )

        # ------------------------------------------
        # plot( total ~ yr, K, subset=which(region=="cfanorth" & vars=="R0.mass" ), type="b" )
        y=list(y=1998, r="cfanorth", v="R0.mass", transgaussian.kriging=T, nruns=1)  # order is important
        # GLS simulation predictions
        PS = kriging.db( DS="UK.conditional.simulation.PS", p=y)
        levelplot( pred.mean~ plon+plat, PS, aspect="iso")

     }




    }

    outdir = file.path( project.datadirectory("bio.snowcrab"), "R", "kriging" )
    dir.create(path=outdir, recursive=T, showWarnings=F)


    if ( DS %in% c( "UK.redo", "UK.point.PS", "UK.conditional.simulation.PS", "UK.conditional.simulation.K", "UK.conditional.simulation.K.complete" ) ) {


      loc.Pnt = file.path( project.datadirectory("bio.snowcrab"), "R", "kriging", "point" )
      loc.sol = file.path( project.datadirectory("bio.snowcrab"), "R", "kriging", "predictions" )
      loc.res = file.path( project.datadirectory("bio.snowcrab"), "R", "kriging", "results" )

      if (p$transgaussian.kriging) {
        loc.Pnt = file.path( loc.Pnt, "trans.gaussian" )
        loc.sol = file.path( loc.sol, "trans.gaussian" )
        loc.res = file.path( loc.res, "trans.gaussian" )
      } else {
        loc.Pnt = file.path( loc.Pnt, "gaussian" )
        loc.sol = file.path( loc.sol, "gaussian" )
        loc.res = file.path( loc.res, "gaussian" )
      }

      dir.create(path=loc.Pnt, recursive=T, showWarnings=F)
      dir.create(path=loc.sol, recursive=T, showWarnings=F)
      dir.create(path=loc.res, recursive=T, showWarnings=F)


      if (DS=="UK.point.PS" ) {
        out = NULL
          v = p$v
          y = p$y
          r = p$r
          fn.Pnt = file.path( loc.Pnt, paste( "PS", v, y, r, "rdata", sep="." ) )
          if ( ! (file.exists( fn.Pnt)) ) return(NULL)
          load(fn.Pnt)
        return ( z )
      }

      if( DS=="UK.conditional.simulation.K.complete") {
          load( p$ofname )
          return( K )
      }

      if (exists( "libs", p)) RLibrary( p$libs )

      if (is.null(ip)) ip = 1:p$nruns

      if (DS %in% c("UK.conditional.simulation.K")  ) {
        out = NULL
        for ( iip in ip ) {
          y = p$runs[iip,"y"]
          r = p$runs[iip,"r"]
          v = p$runs[iip,"v"]
          fn.K = file.path( loc.res, paste( "K", v, y, r, "rdata", sep="." ) )
          if ( ! (file.exists( fn.K)) ) next()
          load(fn.K)

          if (is.null(K$kappa)) K$kappa=NA  ## tmp fix
          if (is.null(K$datestamp)) K$datestamp=NA  ## tmp fix

          out = rbind( out, K )
        }
        K = out
        save( K, file=p$ofname, compress=T )  # glue all the data results together into one file
        return ( K )
      }


      if (DS=="UK.conditional.simulation.PS" ) {
        out = NULL
          v = p$v
          y = p$y
          r = p$r
          fn.PS = file.path( loc.sol, paste( "PS", v, y, r, "rdata", sep="." ) )
          if ( ! (file.exists( fn.PS)) ) return(NULL)
          load(fn.PS)
        return (PS)
      }

      # get data for the whole area
      var.list= c( all.vars( p$kformula), "yr", "plon", "plat")
      K = NULL



      for ( iip in ip ) {

        y = p$runs[iip,"y"]
        r = p$runs[iip,"r"]
        v = p$runs[iip,"v"]

        if ( ! is.null( overwrite.threshold )  ) {
          fn.test= NULL
          if ( "point.kriging" %in% p$kriging.type ) fn.test = file.path( loc.Pnt, paste( "PS", v, y, r, "rdata", sep="." ) )
          if ( "block.conditional.sims" %in% p$kriging.type ) fn.test = file.path( loc.res, paste( "K", v, y, r, "rdata", sep="." ) )
          if ( overwrite.threshold > 0 ) {  # any value <= 0 will force an overwrite
            if ( file.exists( fn.test ) ) {
              load( fn.test)
              current.time = Sys.time()
              file.time = as.POSICct( K$datestamp )
              if (( current.time - file.time ) < overwrite.threshold ) next()
            }
          } else {  # less than 0 -- only overwrite if no results found
            if ( file.exists( fn.test ) ) next()
          }
        }

        S  = snowcrab.db( DS="set.complete")

        if (p$transgaussian.kriging) {
          # S = S[ which( S[,v] > 0) ,]
          S[,v] = variable.recode( S[,v], v, direction="forward", db="snowcrab" )
        }

        S = rename.df(S, v, "kv")
        S = S[, var.list ]

        # jitter is used to prevent singular solutions when stations are too close to each other
        S$plon = jitter(S$plon)
        S$plat = jitter(S$plat)

#        # run a global variogram to obtain constants:: the whole region  (CFAALL and various year combinations are attempted)
#        # the returned semivariance are scaled to unit variance (of the input data)
        vgm.empirical.global = try (model.variogram.running (S, v, y, p, empirical.vgm.only=T), silent=T )
        if (is.null(vgm.empirical.global) | ( "try-error" %in% class(vgm.empirical.global) ) ) {
          print( paste( "Empirical Variogram problem with:", v, y, r ) )
          print (vgm.empirical.global)
          #      next()  # no data for emprical variogram
        }

        iPS = iRegion.PS =iRegion.Sy = z = vgm.e = vgm.m = NULL

        # 1. get data for the year in question
        Sy = S[ which(S$yr==y) ,]
        PS = snowcrab.habitat.db( DS="PS", yrs=y, p=p )


        # 2. habitat area for snow crab
        iHabitat = which( PS$habitat > 0.5 & (PS$habitat.se < PS$habitat ) )  # from a scatterplot, those se > 0.4 seem highly extreme
        # iHabitat = which( PS$habitat > 0.5 & PS$habitat.se < 0.4 )  # from a scatterplot, those se > 0.4 seem highly extreme
        # iHabitat = which( ( PS$habitat - PS$habitat.se*2 ) > 0.5 )

        # 2B ..  # trim PS to make the following less memory intensive as other vars are no longer needed
        PS = PS[ , setdiff(var.list, "kv") ]

        # 3. that which is in the region of interest
        iRegion.Sy = filter.region.polygon(x=Sy[ , c("plon", "plat")], region=r, planar=T)
        iRegion.PS = filter.region.polygon(x=PS[ , c("plon", "plat")], region=r, planar=T)

        totalsurfacearea = length( iHabitat ) * (p$pres*p$pres)  # only habitat space being considered

        # always add areas within a certain distance of a tow, even if they are not in "good habitat"  --- for plotting only
        # keeping all points within X km of a survey station
        distances = rdist( PS[,c("plon", "plat")], Sy[ iRegion.Sy, c("plon", "plat") ])
        dd = distances
        dd[ which(dd < p$threshold.distance) ] =  NA   # mark those within x km of a survey station
        survey.buffer = sort( which( !is.finite( rowSums(dd ) ) ) )

        # plotting surface
        iPlot = union( iHabitat, survey.buffer)
        iPlot = intersect( iRegion.PS, iPlot )

        # prediction surface
        # 5. prediction surface and set-level data for abundance estimation in region of interest
        # including only areas within a given distance from survey stations and east of the SW point of NS
        # remove all areas that are too far from a station p$kriging.extrapolation.limit
        iPS = union( iHabitat, survey.buffer)
        iPS = intersect( iRegion.PS, iPS )
          # before continuing, extract SA estimates
          surfacearea = length(iPS) * (p$pres*p$pres)
        ee = distances
        ee[ which(ee < p$kriging.extrapolation.limit) ] =  NA   # mark those within x km of a survey station
        extrap.buffer = sort( which( !is.finite( rowSums(ee ) ) ) )
        iEastof250 = which( PS$plon > 250 )
        iPS = intersect( iPS, iEastof250 )
        iPS = intersect( iPS, extrap.buffer)

        rm (survey.buffer, dd , distances, iRegion.PS , iHabitat,extrap.buffer, ee  ); gc()

        if ( length(iPS) < 30 ) next() # no available habitat space of significance

        # add a buffer area around the region for input data to reduce edge effects
        i.input.data = iRegion.Sy
        if ( !is.null( p$kriging.distance.buffer ) ) {
          Sdistances =  rdist( Sy[ ,c("plon", "plat")], Sy[ iRegion.Sy, c("plon", "plat") ])
          Sdistances[ which( Sdistances < p$kriging.distance.buffer ) ] =  NA
          S.to.keep = sort( which( !is.finite( rowSums( Sdistances) ) ) )
          i.input.data = union( i.input.data, S.to.keep )
          rm (Sdistances, S.to.keep)
        }
        if (length( i.input.data) < 10 ) next() #no data
        rm( iRegion.Sy ); gc()

        # begin variogram modelling
        if (p$use.local.variograms ) {
          iRegion.S = filter.region.polygon(x=S[, c("plon", "plat")], region=r, planar=T)
          vgm.e0 = NULL
          vgm.e0 = try (model.variogram.running (S=S[ iRegion.S, ], v=v, y=y, p=p, empirical.vgm.only=T), silent=T )
          if (is.null(vgm.e0) | ( "try-error" %in% class(vgm.e0) ) ) {
            vgm.e0 = vgm.empirical.global
            print( paste( "Empirical Variogram problem with:", v, y, r ) )
            print (vgm.e0)
            print( "Using global Variogram instead due to small sample size or unstable variogram" )
          }
        } else {
          vgm.e0 = vgm.empirical.global
        }

        vgm.test = try( gstat.model.variogram ( vgm.e0, vp=get.variogram.params(v, y) ), silent=T)  # the get.vario .. is a blacklist of unusable variables/ years
        if ( "try-error" %in% class(vgm.test) ) {
          print( paste( "Empirical Variogram problem with:", v, y, r ) )
          print (vgm.e0)
          next()  # no data for emprical variogram
        }
        # scale the empirical variogram to the variance in the data in the sub-area and then obtain the variogram model
        var.e = var( Sy[i.input.data,"kv"], na.rm=T )
        vgm.e = vgm.e0  # vgm.e0 is the global variogram for the whole region/year scaled to unit variance
        vgm.e$gamma = vgm.e$gamma * var.e  # rescale variance to that of the specific sub-area/year of interest
        vgm.m = try( gstat.model.variogram ( vgm.e, vp=get.variogram.params(v, y) ), silent=T)  # the get.vario .. is a blacklist of unusable variables/ years

        if ( is.null(vgm.m) | ( "try-error" %in% class(vgm.m) )  ) {
          print( paste( "Modelled variogram problem with:", v, y, r ) )
          print( vgm.e )
          print( vgm.m )
          next()  # no viable solutions
        }

        print( paste( "Kriging with:", v, y, r ) )

        vario.sse = ifelse(is.null(attr(vgm.m, "SSErr")), NA, attr(vgm.m, "SSErr") )
        vario.model = vgm.m$model[2]
        psill = vgm.m$psill[2]
        nugget = vgm.m$psill[1]
        range = vgm.m$range[2]
        kappa = ifelse ( vario.model=="Mat", vgm.m$kappa[2], NA )

        g = gstat(id=v, model=vgm.m, formula=p$kformula, locations=p$klocs, data=Sy[i.input.data,], nmin=1, nmax=p$knmax )
        z = K0 =NULL

        if ( "point.kriging" %in% p$kriging.type ) {
          fn.Pnt = file.path( loc.Pnt, paste( "PS", v, y, r, "rdata", sep="." ) )
          print (fn.Pnt )
         # point kriging .. used for plotting/mapping
          z = try( predict(object=g, newdata=PS[iPlot,], debug=-1), silent=T)
          if ( "try-error" %in% class(z) )  {  # increase the number of stations
             g = gstat(id=v, model=vgm.m, formula=p$kformula, locations=p$klocs, data=Sy[i.input.data,], nmin=1, nmax=p$knmax+20  )
             z = try( predict(object=g, newdata=PS[iPlot,], debug=-1), silent=T)
          }

          if ( "try-error" %in% class(z) )  {  # increase the number of stations
             g = gstat(id=v, model=vgm.m, formula=p$kformula, locations=p$klocs, data=Sy[i.input.data,]  )
             z = try( predict(object=g, newdata=PS[iPlot,], debug=-1), silent=T)
          }
            if ( "try-error" %in% class(z) ) next()

          if (p$transgaussian.kriging)  {
            z[,3] = variable.recode( z[,3], v, direction="backward", db="snowcrab" )
            z[,4] = variable.recode( sqrt( z[,3]), v, direction="backward", db="snowcrab" )
            z[,4] = z[,4] ^ 2
          }
          er = empirical.ranges( db="snowcrab", v )
          z[ which(z[,3]< er[1]) , 3 ] = er[1]  # assume too small to detect .. 0
          z[ which(z[,3]> er[2]) , 3 ] = er[2]  # upper bound
          save (z, file=fn.Pnt, compress=T)
          print( "... Completed successfully, Now plotting ..." )
          map.krige.lattice( M=list(nruns=1, runs=p$runs[iip,], transgaussian.kriging=p$transgaussian.kriging ) )
        }

        if ( "block" %in% p$kriging.type )  {  # not working as expected ...
#         try to use  SpatialPixelsDataFrame
          # block kriging .. old method
#          PS$plon0 = PS$plon
#          PS$plat0 = PS$plat
          centroid = as.data.frame( t( apply( PS, 2, FUN=function(x){ mean (as.numeric(x), na.rm=T)} )))
#          PS$plon = PS$plon0 - centroid$plon
#          PS$plat = PS$plat0 - centroid$plat
          centroid$plon=0
          centroid$plat=0
#          z = try( predict(object=g,  newdata=centroid, block=PS[iPS,c("plon", "plat")], debug=-1), silent=T)
          z = try( predict(object=g,  newdata=centroid, block=PS[iPS,c("plon", "plat")], debug=-1), silent=T)
          if ( "try-error" %in% class(z) )  next()
          datacols = c(3:ncol(z))
          if (p$transgaussian.kriging)  z[,datacols] = variable.recode( z[,datacols], v, direction="backward", db="snowcrab" )
          z.sum = z[,3] * surfacearea
          z.sum.sd = sd( z.sum )
          z.sum.mean = mean ( z.sum )
          ci = quantile( z.sum, probs=c(0.025, 0.5, 0.975), na.rm=T,names=F )
          lbound = ci[1]
          ubound = ci[3]
          median = ci[2]
          K0 = data.frame( yr=y, vars=v, region=r, var=var.e, total=z.sum.mean, median, lbound, ubound,
            totalsurfacearea, surfacearea, vario.model, vario.sse, psill, nugget, range, kappa, datestamp=as.character(Sys.time() ) )
        }

        if (  "block.conditional.sims"  %in% p$kriging.type  ) {
          # block kriging with Gaussian conditional simulation ("trans-Guassian" sensu Cressie)
          # via direct summation of back-transformed data

          fn.PS = file.path( loc.sol, paste( "PS", v, y, r, "rdata", sep="." ) )
          fn.K = file.path( loc.res, paste( "K", v, y, r, "rdata", sep="." ) )

          print(fn.K)

          z = try( predict(object=g, newdata=PS[iPS,], block=c(p$pres,p$pres), nsim=p$n.conditional.sims, debug=-1), silent=T)
          if ( "try-error" %in% class(z) )  { # increase the number of stations
             g = gstat(id=v, model=vgm.m, formula=p$kformula, locations=p$klocs, data=Sy[i.input.data,], nmin=1, nmax=p$knmax+20 )
             z = try( predict(object=g, newdata=PS[iPS,], block=c(p$pres,p$pres), nsim=p$n.conditional.sims, debug=-1), silent=T)
          }
          if ( "try-error" %in% class(z) )  { # increase the number of stations
             g = gstat(id=v, model=vgm.m, formula=p$kformula, locations=p$klocs, data=Sy[i.input.data,], nmin=1, nmax=p$knmax+40 )
             z = try( predict(object=g, newdata=PS[iPS,], block=c(p$pres,p$pres), nsim=p$n.conditional.sims, debug=-1), silent=T)
          }

          if ( "try-error" %in% class(z) )  { # increase the number of stations
             g = gstat(id=v, model=vgm.m, formula=p$kformula, locations=p$klocs, data=Sy[i.input.data,] )
             z = try( predict(object=g, newdata=PS[iPS,], block=c(p$pres,p$pres), nsim=p$n.conditional.sims, debug=-1), silent=T)
          }

          if ( "try-error" %in% class(z) ) next()

          z = as.matrix( z[ , 3:ncol(z) ] )

          # zcd = apply( z,1,mean,na.rm=T ) / apply( z, 1, sd, na.rm=T )
          # ibad = which( zcd > 1 ) # large coefficients of variation are excluded

          if (p$transgaussian.kriging) z = variable.recode( z, v, direction="backward", db="snowcrab" )
          er = empirical.ranges( db="snowcrab", v )

          # z[ ibad,] = 0 # drop them from further consideration

          # remove extremes
          ii = which(z < er[1])
          if (length (ii) > 0 )  z[ ii ] = 0

          jj = which( z > er[2] )
          if (length(jj) > 0 )   z[ jj ] = er[2]

          # remove estimates that are lower than the quantile at which "habitat" is considered viable
          # this is because these low value estimates are highly uncertain but can amount to a large  surface area and so have a large influence
          #ik = which( z < quantile( z[which(z>0)], prob=p$habitat.threshold.quantile ) )
          #if (length (ik) > 0 )  z[ ik ] = 0


          PS = PS[iPS, c("plon", "plat")]
          PS$pred.mean = apply( z, 1, mean, na.rm=T )
          PS$pred.sd = apply( z, 1, sd, na.rm=T )
          save(PS, file=fn.PS, compress=T)

          rm(ii, jj, PS); gc()

          z.sum = apply( z, 2, sum, na.rm=T )
          z.sum.sd = sd( z.sum )
          z.sum.mean = mean ( z.sum )
          ci = quantile( z.sum, probs=c(0.025, 0.5, 0.975), na.rm=T,names=F )
          lbound = ci[1]
          ubound = ci[3]
          median = ci[2]
          K = data.frame( yr=y, vars=v, region=r, var=var.e, total=z.sum.mean, median, lbound, ubound,
            totalsurfacearea, surfacearea, vario.model, vario.sse, psill, nugget, range, kappa, datestamp=as.character( Sys.time() ) )

          save(K, file=fn.K, compress=T)
          print (K)
          print( "... Completed successfully" )
        }

      } # runs (ip)
    } # end UK


  }


