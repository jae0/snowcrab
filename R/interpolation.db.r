
interpolation.db = function( ip=NULL, DS=NULL, p=NULL,
  varnames = c("snowcrab.large.males_abundance", "snowcrab.large.males_presence_absence"), annot.cex=2 ) {

  if (DS %in% c( "fishable.biomass", "fishable.biomass.redo" )) {

    outdir = file.path( project.datadirectory("bio.snowcrab"), "modelled", "biomass" )
    dir.create(path=outdir, recursive=T, showWarnings=F)
    fn = file.path( outdir, "biomass.rdata" )

    if (DS=="fishable.biomass") {
      B = NULL
      if ( file.exists( fn) ) load(fn)
      return (B)
    }

    bm = snowcrab_stmv( p=p, DS="baseline", ret="mean", varnames=varnames )
    bl = snowcrab_stmv( p=p, DS="baseline", ret="lb", varnames=varnames )
    bu = snowcrab_stmv( p=p, DS="baseline", ret="ub", varnames=varnames )

    m = bm[[1]] * bm[[2]]    # biomass (density)
    lb = bl[[1]] * bl[[2]]
    ub = bu[[1]] * bu[[2]]

    h = bm[[2]]  # habitat
    hl = bl[[2]]
    hu = bu[[2]]

    bm=bu=bl = NULL

    # respect the bounds of input data (no extrapolation)
    # set = aegis::survey.db( p=p, DS="filter" ) # mature male > 95 mm
    # qn = quantile( set$totwgt_adjusted, probs=p$stmv_quantile_bounds, na.rm=TRUE )
    # bm[ bm > qn[2] ] = qn[2]  # truncate .. do not extrapolate
    # bm[ bm < qn[1] ] = 0  # these are assumed to be below detection limit
    # bm[ bl < qn[1] ] = 0  # these are assumed to be below detection limit
    bm[ h  < p$habitat.threshold.quantile ] = NA
    bm[ hl < p$habitat.threshold.quantile ] = NA

    if(0) {
      bs = bathymetry.db(p=p, DS="baseline")
      levelplot( m[,16] ~ plon+plat, bs, aspect="iso")
      for (i in 1:16) print(levelplot( m[,i] ~ plon+plat, bs, aspect="iso"))
    }


    # limit range of extrapolation to within a given distance from survey stations .. annual basis
    set = snowcrab.db( DS="set.clean")
    bs = bathymetry.db(p=p, DS="baseline")
    bb = array_map( "xy->1", bs, gridparams=p$gridparams )

    if (0) {
      # annual mask
      for (iy in 1:p$ny ) {
        S = set[ which(set$yr==p$yrs[iy]), c("plon", "plat") ]
        S = S[ !duplicated(S),]
        nn = array_map( "xy->1", S, gridparams=p$gridparams )
        overlap = match(  nn, bb )
        overlap = overlap[ which( is.finite( overlap ))]
        o = bs[overlap,]
        # add corners as buffer
        ot = t(o) # reshape to make addition simpler using R's cycling rules
        o = rbind( t( ot + p$threshold.distance*c( 1, 1) ),
                   t( ot + p$threshold.distance*c(-1,-1) ),
                   t( ot + p$threshold.distance*c( 1,-1) ),
                   t( ot + p$threshold.distance*c(-1, 1) )
        )
        o = o[ !duplicated(o),]
        boundary= non_convex_hull( o, alpha=p$threshold.distance*4, plot=FALSE )
        outside.polygon = which( point.in.polygon( bs[,1], bs[,2], boundary[,1], boundary[,2] ) == 0 )
        m[outside.polygon,iy] = NA
        ub[outside.polygon,iy] = NA
        lb[outside.polygon,iy] = NA
        h[outside.polygon,iy] = NA
        hl[outside.polygon,iy] = NA
        hu[outside.polygon,iy] = NA

      }
    }

    # a static mask
    S = set[ , c("plon", "plat") ]
    S = S[ !duplicated(S),]
    nn = array_map( "xy->1", S, gridparams=p$gridparams )
    overlap = match(  nn, bb )
    overlap = overlap[ which( is.finite( overlap ))]
    o = bs[overlap,]
    # add corners as buffer
    ot = t(o) # reshape to make addition simpler using R's cycling rules
    o = rbind( t( ot + p$threshold.distance*c( 1, 1) ),
               t( ot + p$threshold.distance*c(-1,-1) ),
               t( ot + p$threshold.distance*c( 1,-1) ),
               t( ot + p$threshold.distance*c(-1, 1) )
    )
    o = o[ !duplicated(o),]
    boundary= non_convex_hull( o, alpha=25, plot=FALSE )
    outside.polygon = which( point.in.polygon( bs[,1], bs[,2], boundary[,1], boundary[,2] ) == 0 )
    m[outside.polygon,] = NA
    ub[outside.polygon,] = NA
    lb[outside.polygon,] = NA
    h[outside.polygon,] = NA
    hl[outside.polygon,] = NA
    hu[outside.polygon,] = NA

    B = list( m=m, lb=lb, ub=ub, h=h, hl=hl, hu=hu )
    save( B, file=fn, compress=TRUE )

    return(fn)
  }


  # ------------------


  if (DS %in% c( "fishable.biomass.map" )) {

    projectdir = file.path(p$data_root, "maps", "fishable.biomass", p$spatial.domain )
    dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)

    bs = bathymetry.db(p=p, DS="baseline")
    bm = interpolation.db( p=p, DS="fishable.biomass" )

    fb = bm$m  / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
    fl = bm$lb  / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
    fu = bm$ub  / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
    h  = bm$h

    qs = range(fb[fb>0], na.rm=TRUE)
    datarange = seq( (qs[1]), (qs[2]), length.out=150)
    cols = color.code( "seis", datarange )
    fb[which(!is.finite(fb))] = qs[1]*0.99

    for (iy in 1:p$ny) {
      y = p$yrs[iy]
      xyz = cbind( bs[, c("plon", "plat")], (fb[,iy]) )
      outfn = paste( "prediction.abundance.mean", y, sep=".")
      fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
      png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
      lp = aegis::aegis_map( xyz=xyz, depthcontours=TRUE, pts=NULL, annot=y,
        annot.cex=annot.cex, corners=p$planar.corners, at=datarange,
        col.regions=cols, rez=c(p$pres,p$pres), plotlines="cfa.regions"  )
      print(lp)
      dev.off()
      print (fn)
    }


    qs = range(fl[fl>0], na.rm=TRUE)
    datarange = seq( (qs[1]), (qs[2]), length.out=150)
    cols = color.code( "seis", datarange )
    fl[which(!is.finite(fl))] = qs[1]*0.99

    for (iy in 1:p$ny) {
      y = p$yrs[iy]
      xyz = cbind( bs[, c("plon", "plat")], (fl[,iy]) )
      outfn = paste( "prediction.abundance.lb", y, sep=".")
      fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
      png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
      lp = aegis::aegis_map( xyz=xyz, depthcontours=TRUE, pts=NULL, annot=y,
        annot.cex=annot.cex, corners=p$planar.corners, at=datarange,
        col.regions=cols, rez=c(p$pres,p$pres), plotlines="cfa.regions"  )
      print(lp)
      dev.off()
      print (fn)
    }



    qs = range(fu[fu>0], na.rm=TRUE)
    datarange = seq( (qs[1]), (qs[2]), length.out=150)
    cols = color.code( "seis", datarange )
    fu[which(!is.finite(fu))] = qs[1]*0.99

    for (iy in 1:p$ny) {
      y = p$yrs[iy]
      xyz = cbind( bs[, c("plon", "plat")], (fu[,iy]) )
      outfn = paste( "prediction.abundance.ub", y, sep=".")
      fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
      png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
      lp = aegis::aegis_map( xyz=xyz, depthcontours=TRUE, pts=NULL, annot=y,
        annot.cex=annot.cex, corners=p$planar.corners, at=datarange,
        col.regions=cols, rez=c(p$pres,p$pres), plotlines="cfa.regions"  )
      print(lp)
      dev.off()
      print (fn)
    }

    return(fn)
  }



  # ------------------


  if (DS=="fishable.biomass.timeseries") {
    bm = interpolation.db(p=p, DS="fishable.biomass")

    fb = bm$m  / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
    fl = bm$lb / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
    fu = bm$ub / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
    h  = bm$h

    bs = bathymetry.db( p=p, DS="baseline")

    K = NULL
    nreg = length(p$regions.to.model)
    for (r in 1:nreg ){
      aoi = aegis::polygon_inside(x=bs[ , c("plon", "plat")], region=p$regions.to.model[r], planar=TRUE, proj.type=p$internal.crs )
      aoi = intersect( aoi, which( bs$plon > 250 ) )
      out = matrix( NA, nrow=p$ny, ncol=5)

      for (y in 1:p$ny) {
        iHabitat = which(  {h[,y] >= p$habitat.threshold.quantile  }  ) # any area with biomass > lowest threshold, by definition
        iHabitatRegion = intersect( aoi, iHabitat )
        out[ y, 1] = sum( fb[iHabitatRegion,y] , na.rm=TRUE ) # abundance weighted by Pr
        out[ y, 2] = sum( fl[iHabitatRegion,y] , na.rm=TRUE )
        out[ y, 3] = sum( fu[iHabitatRegion,y] , na.rm=TRUE )
        out[ y, 4] = sum( h[iHabitatRegion,y] ) * (p$pres*p$pres)
        out[ y, 5] = length( iHabitatRegion ) * (p$pres*p$pres)
      }

      ok = as.data.frame( out )
      names( ok) = c("total", "total.lb", "total.ub", "sa.region", "sa.crude")

      ok$log.total = log(ok$total)
      ok$log.total.lb = log( ok$total.lb) # as above
      ok$log.total.ub = log( ok$total.ub) # as above
      ok$region = p$regions.to.model[r]
      ok$yr = p$yrs

      K = rbind(K, ok)
    }

    K$density = K$total / K$sa.region
    K$density.crude = K$total / K$sa.crude

    return( K )

    if (0){
      str(K)
      table.view( K )
      plot( total ~ yr, K[K$region=="cfanorth", ], type="b")
      plot( total ~ yr, K[K$region=="cfasouth", ], type="b")
      plot( total ~ yr, K[K$region=="cfa4x", ], type="b")

      plot( density ~ yr, K[K$region=="cfanorth", ], type="b")
      plot( density ~ yr, K[K$region=="cfasouth", ], type="b")
      plot( density ~ yr, K[K$region=="cfa4x", ], type="b")

      plot( density.crude ~ yr, K[K$region=="cfanorth", ], type="b")
      plot( density.crude ~ yr, K[K$region=="cfasouth", ], type="b")
      plot( density.crude ~ yr, K[K$region=="cfa4x", ], type="b")

    }

  }


  if ( DS %in% c( "interpolation.simulation" ) ) {
    message(" simulation-based results are not ready at present")
    message(" defaulting to simple estimates based upon asymptotic assumptions" )
    message(" only R0.mass is supported for now" )

    out = NULL
    if ( p$vars.to.model == "R0.mass" ) out = interpolation.db( p=p, DS="fishable.biomass.timeseries" )

    return(out)
  }


# ---------


  if (DS =="habitat.temperatures") {

    bm = interpolation.db( p=p, DS="fishable.biomass" )
    ps = snowcrab_stmv(p=p, DS="output_data" )
    bs = bathymetry.db( p=p, DS="baseline")

    temp = ps$t

    K = NULL
    nreg = length(p$regions.to.model)
    for (r in 1:nreg ){
      aoi = aegis::polygon_inside(x=bs[ , c("plon", "plat")], region=p$regions.to.model[r], planar=TRUE, proj.type=p$internal.crs )
      aoi = intersect( aoi, which( bs$plon > 250 ) )
      out = matrix( NA, nrow=p$ny, ncol=2)

      for (y in 1:p$ny) {
        iHabitat = which( bm$h[,y] >= p$habitat.threshold.quantile ) # any area with bm > lowest threshold
        iHabitatRegion = intersect( aoi, iHabitat )
        out[ y, 1] = mean( temp[iHabitatRegion,y] , na.rm=TRUE ) # temperature weighted by Pr
        out[ y, 2] = sd( temp[iHabitatRegion,y] , na.rm=TRUE ) # temperature weighted by Pr
      }

      ok = as.data.frame( out )
      names( ok) = c("temperature", "temperature.sd")
      ok$region = p$regions.to.model[r]
      ok$yr = p$yrs
      ok$lbound = ok$temperature - ok$temperature.sd*1.96 # normal assumption
      ok$ubound = ok$temperature + ok$temperature.sd*1.96
      K = rbind(K, ok)
    }

    return( K )

  }


  return ("Completed" )
}
