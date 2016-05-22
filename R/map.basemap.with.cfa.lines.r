
  map.basemap.with.cfa.lines = function( p, conversions=c("ps2png") ) {
    # map blank with cfa lines
    p$tension = "-T.4"  # 0.35+ for steep; 0.25 for smooth
    p$maskres = "-S16k"
    p$interpres = "-nb"
    set = snowcrab.db( DS="set.clean")
    set$dummy = 1
    set$sa = 1
    u = set[, c("lon", "lat", "dummy", "sa")]
    u = u[is.finite(u[,"lon"]*u[,"lat"] ),]
    dir.create ( p$outdir, recursive=T, showWarnings=F )
    toplot = u[1,]
		tmpdir  = tempdir()
    with( p, {
    # tmp local GMT files
      grid  = file.path(tmpdir, make.random.string("gmt.grid"))
      bin   = file.path(tmpdir, make.random.string(".gmt.bin"))  # binary data
      clip  = file.path(tmpdir, make.random.string(".gmt.clip"))
      indat = file.path(tmpdir, make.random.string(".gmt.indat"))
      blocked = file.path(tmpdir, make.random.string(".gmt.blocked"))
      annotations = file.path(tmpdir, make.random.string(".gmt.annotations"))
      outfile = paste( file.path(tmpdir, make.random.string(".gmt.outfile")), ".ps", sep="") # final postscript file
      write.table( toplot, file=indat, quote=F, col.names=F, row.names=F)
    # define some local parameters
      append = "-O -K"
      gmt.prep (fnt=fnt, psres=psresolution)
    # begin mapping
      cmd( "psbasemap", region, gmtproj, annot, "-K >", outfile )
    # various overlays
      cmd( "cat", basemap, ">>", outfile ) # bathymetry
      if (!is.null(overlay)) {
      if ( length(overlay)>0 & overlay != "" ) {
        for (o in overlay) {
          cmd( "psxy", find.ecomod.gis(o), region, gmtproj, append, polygon.options, ">>", outfile ) # polygons
      }}}
      cmd( "pscoast", region, gmtproj, coast.options, "-O >>", outfile ) # coastline
      remove.files ( c( grid, bin, clip, indat, blocked, annotations ))
      if (!is.finite(conversions)) { #assume postscript not needed
        graphics.convert ( outfile, p$outfile.basename, conversions )
      } else { # assume just the postscript file is required
        file.copy(from=outfile, to=paste(p$outfile.basename,"ps",sep="."), overwrite=T)
      }
    }) # end with
    gmt.cleanup ()
    return("Done") 
  }

