
  figure.timeseries.interpolated = function( p, outdir ) {
    
    p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
    K = interpolation.db( DS="interpolation.simulation", p=p  ) 
    
    for (v in p$vars.to.model) {
    for (r in p$regions.to.model) {
      print (paste(v,r))
      i = which(K$vars==v & K$region==r )
      if (length(i) == 0) next()
      XX=K[,]
      xx = XX$yr
      yy = XX$total
      lbound = XX$lbound
      ubound = XX$ubound
      surfacearea = XX$surfacearea

      dir.create( outdir, recursive=T, showWarnings=F  )
      ts.plotandsave(xx, yy, lbound, ubound, surfacearea, action="save", title="", filename=paste(v,r,sep="."), outdir=outdir )
    }}
    return("Done") 
  }



