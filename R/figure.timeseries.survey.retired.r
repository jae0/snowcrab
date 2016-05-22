  figure.timeseries.survey.retired = function(p, outdir=NULL, areas=NULL, from.file=F) {
     
    if (is.null(areas)) areas = c( "cfa20", "cfa21", "cfa22", "cfa23", "cfa24", "cfa4x", 
      "cfa23slope", "cfa24slope", "cfaslope", "cfanorth", "cfasouth", "cfaall" )
    variables =  variable.list.expand("all.data")
    
    if (from.file) {
      tsdata =  get.time.series ( from.file=T )  # this returns 1.96SE as "se"
    } else {
      tsdata =  get.time.series ( x=snowcrab.db( DS="set.logbook" ), regions=areas, vars=variables, from.file=F, trim=0 )
    }

    plot.timeseries(x=tsdata, vars=variables, regions=areas, outdir=outdir, backtransform=F)
#    unloggedvars =  c("sexratio.all", "sexratio.mat", "sexratio.imm", "z", "t", "julian")

#    v1 = setdiff(variables, unloggedvars)
#    if (length(v1) > 0) plot.timeseries(x=tsdata, vars=v1, regions=areas, outdir=outdir, backtransform=T)
    
#    v2 = intersect(variables, unloggedvars )
#    if (length(v2) > 0) plot.timeseries(x=tsdata, vars=v2, regions=areas, outdir=outdir, backtransform=F)
    
    return("Done")
  }


