 
  create.isobath.data = function(p, depth=1000, outdir=polygonsdir, outfile="isobath1000m.dat") {
    d = gmt.isobath (p, depth=depth)
    d = d[,c(1,2)]
    plot(d$V1, d$V2, type="b")
    if (depth==1000) {
      pt1 = cbind(-56, 45)
      pt2 = cbind(-68, 42)
      pt3 = cbind(-66, 44)
      pt4 = cbind(-61, 47)
      pt5 = cbind(-60, 48)
      d = rbind(pt1, d, pt2, pt3, pt4, pt5, pt1)
    }
    write.table(d, file.path( outdir, outfile, col.names=F, row.names=F, quote=F))
    return(d)
  }

