
      map.surface = function( ip=NULL, plotvar, pyears, vclass, log.transf=F, ... ) {

      # ip is the first parameter passed in the parallel mode

        if (exists( "libs", p)) RLibrary( p$libs )
        if (is.null(ip)) ip = 1:length(pyears)

        outdir = file.path( project.datadirectory("bio.snowcrab"), "R", "habitat", "maps" )
        dir.create(path=outdir, recursive=T, showWarnings=F)

        if ( plotvar %in% c("habitat", "habitat.sim") ) {
          dr = c( 0, 1)
        }
        if ( plotvar %in% c("abundance", "abundance.sim") ) {
          xx = snowcrab.db("set.biologicals")[, "R0.mass" ]
          dr = range( xx, na.rm=T )
          if (log.transf) {
            er = range( xx[ which(xx >0)])
            dr = log10(er)
          }
        }

        datarange = seq( dr[1] , dr[2], length.out=100 )
        cols = color.code( "seis", datarange )

        for (iy in ip) {
          y = pyears[iy]
          fn =  paste( plotvar, y, vclass, sep=".")
          print ( file.path(outdir,fn) )
          PS =   predictive.simulations( DS="PS", vclass=vclass, pyears=y )

          if (log.transf) {
            PS[, plotvar ] = log10( PS[, plotvar ] )
          }
          map( PS[, c("plon", "plat", plotvar)], xyz.coords="planar", cfa.regions=T, depthcontours=T, annot=y,
            fn=fn, loc=outdir, at=datarange, col.regions=cols, corners=planar.corners)  # planar.corners defined in load.environment.r
        }

        return( "Complete" )
      }


