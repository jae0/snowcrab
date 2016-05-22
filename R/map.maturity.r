
  map.maturity = function(p, outdir ) {
    
    load(file.path( project.datadirectory("snowcrab"), "R", "maturity.rdata"))
    
    p$tension = "-T.4"  # 0.35+ for steep; 0.25 for smooth
    p$maskres = "-S16k"
    p$interpres = "-nb"

    for (sex in unique(maturity$sex)) {
      x = maturity[ which(maturity$sex==sex) ,]
      x$sa =  1

      variables = c("cw50")
      x = x[ is.finite(rowSums(x[, c("yr", "lon", "lat")])), ]
      basedir = "maturity"
      outdir = file.path( basedir, paste("sex",sex,sep=""), p$spatial.domain )
      gmt.map.variables( x, p, variables, plottimes=p$plottimes, outdir, p$conversions, init.files=p$init.files )
    }
  }


