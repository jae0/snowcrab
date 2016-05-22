 plot.fisheries.ts.retired = function(x, variable, directory) {
    datarange = range(x$yr, na.rm=T) + c(-1, 1)
    regions  = sort(unique(x$cfa))
    for (reg in regions) {
      u = x[ which(x$cfa==reg), ]
      u = u[ which(is.finite(u[,variable])),]
      u = u[order(u$yr),]
      #if (dim(u)[1]>3) {
      #outfile = paste( variable, reg, sep="" )
      plot( u$yr, u[,variable], type="b", axes=F, xlab="Years",
              ylab=paste(variable, reg, sep="."), xlim=datarange, col="orange", lty="solid", lwd=5 )
        points(x=u$yr, y= u[,variable], pch=10)
        axis( 1 )
        axis( 2 )

      Pr(dev="png", dname=directory, fname=outfile, width=8, height=4)

      print (outfile)
      }
   }
    return(NULL)

  }

