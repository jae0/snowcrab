
plot.projections = function(p, pj, reg, outdir, reg.label) {
  default.options = options()
  setup.lattice.options()
  xlim=range(pj$yr); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5
  ylim=range(pj$ub); ylim[1]=0; ylim[2]=  ylim[2]+ 0.5
  nex = length( unique ( pj$ER ) )
  uu = xyplot( fbpre~yr|ER, data=pj, ub=pj$ub, lb=pj$lb, reg.label=reg.label, 
            layout=c(1,nex), xlim=xlim, ylim=ylim,
            main=paste("Projections of fishable biomass relative to", p$start.projection.year, "--", reg ),
            ylab = "Fishable biomass (t; pre-fishery)", xlab= "Year",
#            main = reg.label, 
#            xlab="Year", ylab="Relative prefishery fishable biomass (proportion)",
            panel = function(x, y, subscripts, ub, lb, ...) {
#           panel.abline(h=y[1], col="gray75", ...)
            panel.abline(h=1, col="gray75", ...)
            larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length=0.05 )
            panel.xyplot(x, y, type="b", lty="dashed", lwd=1.4, pch=19, col="black", ...)
         }
    )
    plot(uu)
   Pr( dev="png", dname=outdir, fname=paste("projections", reg, sep="."), width=6, height=8)
   options(default.options)
}


