
figure.timeseries.bycatch.jonahcrab = function(outdir, all.areas = T) {
  #browser()
  
  if (all.areas) {
    areas = c("cfa4x", "cfasouth", "cfanorth")
    regions = c("4X", "S-ENS", "N-ENS")
  } else {
    areas = c("cfasouth", "cfanorth")
    regions = c("S-ENS", "N-ENS")
  }
  
  n.regions = length(regions)
  n.areas = length(areas)
  
  v = "ms.mass.2511"
  #MG make a quick recode for these so it can loop through and print the species names based on recode info.
  #cod, halibut, thornyskate, wolfish, lessertoadcrab, jonahcrab, smoothskate, winterskate, northernshrimp,
  #species = c(10, 30, 201, 50, 2521, 2511, 202, 204, 2211)
  
  td =  get.time.series (
    from.file = T, outfile = file.path(project.datadirectory("snowcrab"), "R", "tsbycatch.rdata")
  )
  td = td[which(td$variable == v) ,]
  td = td[order(td$region, td$year) ,]
  td = td[which(td$year > 2003),]
  td$region = factor(td$region, levels = areas, labels = regions)
  #   td[ which(td$region=="4X" & td$year < 2004), c("mean", "se", "ub", "lb", "n")] = NA
  
  ylim = 'NULL'
  ylim[2] = max(td$ub) * 0.8
  ylim[1] = round(min(td$lb), 0) - (max(td$ub) * 0.05)
  #ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.1*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
  xlim = range(td$year); xlim[1] = xlim[1]; xlim[2] = xlim[2]
  
  xlabels = seq(min(xlim), max(xlim), 1)
  ylabels = round(seq(0, round(as.numeric(ylim[2]),-1), length.out = 8),-1)
  
  dir.create(outdir, recursive = T, showWarnings = F)
  fn = file.path(outdir, paste(v, "pdf",  sep = "."))
  cex.main = 1.4
  cex.lab = 1
  cex.axis = 0.2
  pdf(
    file = fn, width = 5, height = 6.5, bg = 'white'
  )
  
  #Cairo( file=fn, type="png", bg="white",  units="in", width=5, height=6.5, dpi=350 )
  setup.lattice.options()
  pl = xyplot(
    mean ~ year | region, data = td, ub = td$ub, lb = td$lb,
    layout = c(1,n.regions),
    par.strip.text = list(
      plot.symbol = list(
        col = 'black', fill = 'darkgrey', cex = 0.75, pch = 21
      ),
      axis.text = list(cex = 0.7),
      par.main.text = list(cex = 1),
      layout.heights = list(
        strip = 1, panel = 1, main = 0.5
      ),
      strip.background = list(col = 'lightgrey')
    ),
    #xlim=xlim,
    ylim = (c(as.numeric(ylim[1]), as.numeric(ylim[2]))),
    scales = list(
      y = list(at = ylabels, labels = ylabels, cex = 0.65), x = list(
        at = xlabels, labels = xlabels, rot = 50, cex = 0.65
      )
    ),
    main = "Jonah Crab Biomass", xlab = list("Year", cex = 1), ylab =
      list("Geometric mean kg / km^2", cex = 1),
    #cex.lab=cex.lab,
    cex.axis = cex.axis,
    cex.main = cex.main,
    panel = function(x, y, subscripts, ub, lb, ...) {
      larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length =
                0.05)
      panel.abline(
        h = median(y), lty = 1, lwd = 1, col = "grey27", ...
      )
      panel.xyplot(
        x, y, type = "b", lty = 1, lwd = 1.5, pch = 21, fill = 'darkgrey', col =
          "black", ...
      )
      
    }
  )
  
  print(pl)
  dev.off()
  #cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
  return("Done")
}
