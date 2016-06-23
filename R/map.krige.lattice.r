  map.krige.lattice = function( ip=NULL, M=NULL, log.transf=T,
    loc = file.path( project.datadirectory("bio.snowcrab"), "R", "kriging", "maps" ) ,
    dat = snowcrab.db("set.complete")
    ) {

    if (M$transgaussian.kriging) {
      loc = file.path( loc, "trans.gaussian" )
    } else {
      loc = file.path( loc, "gaussian" )
    }

    dir.create(path=loc, recursive=T, showWarnings=F)

    # ip is the first parameter passed in the parallel mode
    if (exists( "libs", p)) RLibrary( p$libs )
    if (is.null(ip)) ip = 1:M$nruns

    for (iip in ip) {
        y = M$runs[iip,"y"]
        r = M$runs[iip,"r"]
        v = M$runs[iip,"v"]

      k.pname = paste(v, "pred", sep=".")
      k.vname = paste(v, "var",  sep=".")

      PS = kriging.db( DS="UK.point.PS", p=list(v=v, y=y, r=r, transgaussian.kriging=M$transgaussian.kriging) )
      if (is.null(PS)) next()

      xx = dat[, v]
      er = range( xx, na.rm=T )

      if (log.transf) {
        er = range( xx[ which(xx >0)])
        er = log10(er)
        PS[, k.pname] = log10( PS[,k.pname] )
        PS[ which(PS[, k.pname] < er[1]) , k.pname] = er[1]
        PS[, k.vname] = log10( sqrt(PS[,k.vname])  )
      }

      # mean estimate
      datacols = c("plon", "plat", k.pname)
      datarange = seq( er[1], er[2], length.out=150)
      cols = color.code( "seis",datarange )
      outfn = paste( "prediction.mean", v, y, r, sep=".")
      annot = paste( v, y )
      map( xyz=PS[,datacols], cfa.regions=T, depthcontours=T, pts=NULL, annot=annot,
          annot.cex=M$annot.cex, corners=planar.corners ,
        fn=outfn, loc=loc, at=datarange , col.regions=cols )

      # SD estimate
      datacols = c("plon", "plat", k.vname)
      datarange = seq( er[1]/2, er[2]*2, length.out=150)
      cols = color.code( "seis",datarange )
      outfn = paste( "prediction.sd", v, y, r, sep=".")
      annot = paste( v, y)
      map( xyz=PS[,datacols], cfa.regions=T, depthcontours=T, pts=NULL, annot=annot,
          annot.cex=M$annot.cex, corners=planar.corners ,
        fn=outfn, loc=loc, at=datarange , col.regions=cols )

    }
    return ("Complete")
  }


