
  ts.landings.and.tacs = function( x, y, z, outdir, outfile, title) {
    # one time thing .. silly but faster this way
      plot(x, y, type="n", axes=F, xlab="Years", ylab=title, cex=2 )
      lines( x, y, col="orange", lty="solid", lwd=4 )
      points(x, y, pch=10)

#    lines( x, z, col="blue", lty="solid", lwd=4 )
#      points(x, z, pch=11)

      axis( 1 )
      axis( 2 )

      Pr(dev="png", dname=outdir, fname=outfile, width=8, height=4)

  }

