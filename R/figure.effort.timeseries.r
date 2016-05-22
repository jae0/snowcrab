figure.effort.timeseries = function(yearmax, outdir = NULL, outfile = NULL, outfile2 =
                                      NULL, type = "line") {
  regions = c("cfanorth", "cfasouth", "cfa4x")
  e = NULL
  for (r in regions) {
    res = get.fishery.stats.by.region(Reg = r)
    e = cbind(e, res$effort)
  }
  
  e = e / 1000
  
  e = as.data.frame(e)
  colnames(e) = regions
  rownames(e) = res$yr
  
  e = e[which(as.numeric(rownames(e)) <= yearmax),]
  uyrs = as.numeric(rownames(e))
  
  
  dir.create(outdir, recursive = T, showWarnings = F)
  fn = file.path(outdir, paste(outfile,"pdf",sep = "."))
  # Cairo( file=fn, type="png", bg="white",, pointsize=30, units="in", width=6, height=4, dpi=300  )
  fn2 = file.path(outdir, paste(outfile2,"pdf",sep = "."))
  
  #png( file=fn,units='in', width=7,height=7,pointsize=10, res=350,type='cairo')
  pdf(file = fn, width = 7, height = 7, bg = 'white')
  
  if (type == "bar") {
    e[is.na(e)] = 0
    formed.data = t(as.matrix(e))
    barplot(
      formed.data, space = 0, xlab = "Year", ylab = "Effort (1000 trap hauls)", col =
        cols
    )
    legend(
      x = 1, y = 130, c("N-ENS", "S-ENS", "4X"), fill = cols[reverse], bty = "n"
    )
  }
  if (type == "line") {
    pts = c(19, 22, 24)
    lns = c(1, 1, 1)
    cols = c("grey10", "grey10",  "grey20")
    e[which(e == 0)] = NA
    yrange = range (e, na.rm = T)
    yrange[1] = 0
    xrange = range(uyrs)
    xrange[1] = xrange[1]
    xrange[2] = xrange[2]
    xlabels = seq(xrange[1] + 1, xrange[2], 2)
    
    
    m = 1; plot(
      uyrs, e[,m],  type = "b", ylab = "Effort (1000 trap hauls)", xlab = "Year", col =
        cols[m], lwd = 3, lty = lns[m], pch = pts[m], xaxt = "n", xlim = xrange, ylim =
        yrange
    )
    m = 2; points(
      uyrs, e[,m], type = "b", col = cols[m], lwd = 3, lty = lns[m], pch = pts[m]
    )
    m = 3; points(
      uyrs, e[,m], type = "b", col = cols[m], lwd = 3, lty = lns[m], pch = pts[m]
    )
    axis(1, at = xlabels, labels = FALSE)
    text(
      x = xlabels + 1, y = par('usr')[3], labels = xlabels, srt = 45, adj = c(1.5,1), xpd =
        TRUE
    )
    axis(2)
    legend(
      x = 1980, y = 100, c("N-ENS", "S-ENS", "4X"), bty = "n", lty = lns, lwd =
        2, pch = pts, col = cols, cex = 1.2
    )
    dev.off()
    
   # png(file = fn2 ,units = 'in', width = 7, height = 7,pointsize = 10, res = 350,type ='cairo')
   pdf(file=fn2, width=7, height=7, bg='white')
    
    sm = e[, c(1, 3)]
    pts = c(19, 24)
    lns = c(1, 1)
    cols = c("grey10", "grey20")
    yrange = range (sm, na.rm = T)
    yrange[1] = 0
    xrange = range(uyrs)
    xrange[1] = xrange[1]
    xrange[2] = xrange[2]
    xlabels = seq(xrange[1] + 1, xrange[2], 2)
    
    
    m = 1; plot(
      uyrs, sm[,m],  type = "b", ylab = "Effort (1000 trap hauls)", xlab = "Year", col =
        cols[m], lwd = 4, lty = lns[m], pch = pts[m], xaxt = "n", xlim = xrange, ylim =
        yrange
    )
    m = 2; points(
      uyrs, sm[,m], type = "b", col = cols[m], lwd = 3, lty = lns[m], pch = pts[m]
    )
    axis(1, at = xlabels, labels = FALSE)
    text(
      x = xlabels + 1, y = par('usr')[3], labels = xlabels, srt = 45, adj = c(1.5,1), xpd =
        TRUE
    )
    axis(2)
    legend(
      x = 1985, y = 45, c("N-ENS", "4X"), bty = "n", lty = lns, lwd = 2, pch =
        pts, col = cols, cex = 1.2
    )
  }
  dev.off()
 # cmd("convert -trim -frame 10x10 -mattecolor white ", fn, fn)
  #table.view( e)
  return(fn)
}
