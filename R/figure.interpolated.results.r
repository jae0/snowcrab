
  figure.interpolated.results = function( p, outdir, all.areas=T, alt.zero.y=F ) {
    
    varstoplot = matrix(ncol=3, byrow=T, data=
      c( "R0.mass", "Fishable biomass", "B",
       "R0.no", "Number of fishable crab", "N",
       "R1.no", "Number of recruits", "N",

       "totno.female.imm",  "Number of immature females", "N",
       "totno.female.mat",  "Number of mature females", "N",
       "totno.female.berried",  "Number of berried females", "N",
       "totno.female.primiparous", "Number of primiparous females", "N",
       "totno.female.multiparous",  "Number of multiparous females", "N",

       "mi6.no", "Number of immature - instar 6", "N",
       "mi7.no", "Number of immature - instar 7", "N",
       "mi8.no", "Number of immature - instar 8", "N",
       "mi9.no", "Number of immature - instar 9", "N",
       "mi10.no", "Number of immature - instar 10", "N",
       "mi11.no", "Number of immature - instar 11", "N",
       "mi12.no", "Number of immature - instar 12", "N",

       "mi9.skip.moulter.no", "Number of immature skip moulter - instar 9", "N",
       "mi10.skip.moulter.no", "Number of immature skip moulter - instar 10", "N",
       "mi11.skip.moulter.no", "Number of immature skip moulter - instar 11", "N",
       "mi12.skip.moulter.no", "Number of immature skip moulter - instar 12", "N",

       "ma9.CC1to2.no", "Number of mature - instar 9 CC1 & 2", "N",
       "ma10.CC1to2.no", "Number of mature - instar 10 CC1 & 2", "N",
       "ma11.CC1to2.no", "Number of mature - instar 11 CC1 & 2", "N",
       "ma12.CC1to2.no", "Number of mature - instar 12 CC1 & 2", "N",
       "ma13.CC1to2.no", "Number of mature - instar 13 CC1 & 2", "N",

       "ma9.CC3to4.no", "Number of mature - instar 9 CC3 & 4", "N",
       "ma10.CC3to4.no", "Number of mature - instar 10 CC3 & 4", "N",
       "ma11.CC3to4.no", "Number of mature - instar 11 CC3 & 4", "N",
       "ma12.CC3to4.no", "Number of mature - instar 12 CC3 & 4", "N",
       "ma13.CC3to4.no", "Number of mature - instar 13 CC3 & 4", "N",

       "ma9.CC5.no", "Number of mature - instar 9 CC5", "N",
       "ma10.CC5.no", "Number of mature - instar 10 CC5", "N",
       "ma11.CC5.no", "Number of mature - instar 11 CC5", "N",
       "ma12.CC5.no", "Number of mature - instar 12 CC5", "N",
       "ma13.CC5.no", "Number of mature - instar 13 CC5", "N"

    ) )


    eps = 1e-6

    p$vars.to.model = varstoplot[,1]
    p$years.to.model = 1998:p$current.assessment.year
    p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
    K = interpolation.db( DS="interpolation.simulation", p=p )

    current.year = p$current.assessment.year
    K = K[ -which( as.character(K$region)=="cfa4x" & K$yr <= 2002 ), ]
    K = K[ which( K$yr >1997 ), ]


    if (all.areas) {
      areas = c("cfa4x", "cfasouth", "cfanorth" )
      regions = c("4X", "S-ENS", "N-ENS")
    } else {
      areas = c("cfasouth", "cfanorth" )
      regions = c("S-ENS", "N-ENS")
    }

    n.regions = length(regions)
    K$region = factor(as.character(K$region), levels=areas, labels=regions)

    zero.value = K[1,]
    zero.value$yr = current.year
    zero.value$total = 0
    zero.value$lbound = 0
    zero.value$ubound = 0


    xlim=range(K$yr, na.rm=T)
    xlim[1]=xlim[1]-0.5
    xlim[2]=xlim[2]+0.5
  
  
    nvarstoplot = dim(varstoplot)[1]

    for (i in c(1:nvarstoplot)) {

      v = varstoplot[i,1]
      tt = varstoplot[i,2]

      td = K[ which(as.character(K$vars)==v) ,]

      if (varstoplot[i,3] == "B") {
        yy = expression(paste( "Biomass (X ", 10^3, " t)"))
        convert = 10^3 # from t to kt 
        eps = 1e-6
      } else if (varstoplot[i,3] == "N") {
        yy = expression(paste( "Number (X ", 10^6, " )"))
        convert = 1  # convert to
        eps = 1e-6
      }

      if (v=="R0.mass") {
        ii = which(td$yr==1998 & td$region=="S-ENS")
        if ( length(ii)>0 ) td[ ii, "total" ] = NA  # strange data
      }
      if (v=="R1.no") {
  #      mm = zero.value
  #      mm$region = "N-ENS"
  #      td = rbind(td, mm)
      }

      varstocheck = c("total", "ubound", "lbound")
      for (vs in varstocheck) {
        td[,vs] = td[,vs] / convert
        kk = which(td[,vs] <= eps)
        if ( length(kk)>0 ) td[kk,vs] = 0
      }

      td = td[is.finite(td$total) ,]
      td = td[order(td$region,td$yr),]

      # last check to make sure data exists for all regions for the plots
      limits = list()
      for (rr in 1:n.regions) {
        nr = NULL
        pp = which( td$region==regions[rr] )
        nr = length( pp )
        if (nr ==0) {
          mm = zero.value
          mm$region = regions[rr]
          td = rbind( td, mm)
          qq = 1
        } else {
          qq = max(td$ubound[pp], na.rm=T) 
        }
        limits[[rr]] = c(0, qq )
      }
     
      dir.create( outdir, recursive=T, showWarnings=F  )
      fn = file.path( outdir, paste(v, "png", sep=".") )
      print (fn) 

      Cairo( file=fn, type="png", bg="white", units="in", width=5, height=6, dpi=300)
      setup.lattice.options()
     
      if (alt.zero.y) {

        pl = xyplot( total~yr|region, data=td, ub=td$ubound, lb=td$lbound,
          layout=c(1,n.regions), xlim=xlim, scales = list(y = list(relation="free", limits=limits  )),
          par.strip.text = list(cex=3),
          par.settings=list(  
            axis.text=list(cex=2.5), 
            par.main.text = list(cex=4),
            layout.heights=list(strip=0.2, panel=1, main=0.5 ) 
          ),
              main=list(label=tt), xlab=list(label="Survey year", cex=3.5), ylab=list(label=yy, cex=3.5), 
              panel = function(x, y, subscripts, ub, lb, ...) {
              larrows(x, lb[subscripts],
                      x, ub[subscripts],
                     angle = 90, code = 3, length=0.02, lwd=6)
              panel.xyplot(x, y, type="b", lty=1, lwd=6, pch=19, col="black", cex=1.2, ...)
              panel.abline(v=2001.5, col="gray", lty=5, lwd=4, ...)
              panel.abline(h=0, col="gray", lty=1, lwd=4, ...)
              }
            )
      } else {
        pl = xyplot( total~yr|region, data=td, ub=td$ubound, lb=td$lbound,
          layout=c(1,n.regions), 
          par.strip.text = list(cex=1),
          par.settings=list(  
            axis.text=list(cex=2.5), 
            par.main.text = list(cex=4),
            layout.heights=list(strip=0.2, panel=1, main=0.5 ) 
          ),
          xlim=xlim, scales = list(y = list(relation="free") ),
              main=tt, xlab="Survey year", ylab=list(yy, cex=3), 
              panel = function(x, y, subscripts, ub, lb, ...) {
              larrows(x, lb[subscripts],
                      x, ub[subscripts],
                     angle = 90, code = 3, length=0.05)
              panel.xyplot(x, y, type="b", lty=1, lwd=6, pch=30, col="black", ...)
              panel.abline(v=2001.5, col="black", lty=3, lwd=5, ...)
              panel.abline(h=0, col="black", lty=1, lwd=5, ...)
              }
            )
      } 
     
     print(pl)
     dev.off()
     cmd( "convert -trim -frame 10x10 -mattecolor white ", fn, fn )
     
    }
    return("Done")
  }  


