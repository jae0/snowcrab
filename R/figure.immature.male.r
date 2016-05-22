
  figure.immature.male = function( p, outdir, all.areas=T ) {

    eps = 1e-6
    classes = paste( "mi", c(12:8), ".no", sep="")
    classnames = paste( "Instar", c(12:8) )

    p$vars.to.model = classes
    p$years.to.model = 1998:p$current.assessment.year
    p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
    K = interpolation.db( DS="interpolation.simulation", p=p )

    current.year = p$current.assessment.year
    K = K[ -which( as.character(K$region)=="cfa4x" & K$yr <= 2002 ), ]
    K = K[ -which( K$yr <= 1998 ), ]

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

    cex.main = 1.4
    cex.lab = 1.3
    cex.axis = 1.3

    xlim=range(K$yr, na.rm=T)
    xlim[1]=xlim[1]-0.5
    xlim[2]=xlim[2]+0.5

    n.panels = length(classes)

    for (i in c(1:n.regions)) {

      ll = which( as.character(K$vars) %in% classes  & as.character(K$region) == regions[i] ) 
      if (length(ll) == 0 ) next()
      td = K[ ll ,]
      td$vars = factor( as.character(td$vars), levels=classes, labels=classnames )

      convert = 1  # convert  millions to ...
      eps = 1e-6

      varstocheck = c("total", "ubound", "lbound")
      for (vs in varstocheck) {
        td[,vs] = td[,vs] / convert
        kk = which( td[,vs] <= eps )
        td[kk,vs] = 0
      }

      td = td[is.finite(td$total) ,]
      td = td[order(td$vars,td$yr),]

     # last check to make sure data exists for all regions for the plots
      for (rr in 1:length(classes)) {
        nr=NULL
        nr = length( which(td$vars==classnames[rr] ))
        if (nr ==0) {
          mm = zero.value
          mm$vars = classnames[rr]
          td = rbind( td, mm)
        }
      }

     dir.create( outdir, recursive=T, showWarnings=F  )
     fn = file.path( outdir, paste( "imm",areas[i], sep="." ) )
     Cairo( file=fn, type="pdf", bg="white", units="in", width=6, height=4)
      setup.lattice.options()
      pl = xyplot( total~yr|vars, data=td, ub=td$ubound, lb=td$lbound,
        layout=c(1,n.panels), xlim=xlim, scales = list(y = "free"),
            main=paste("Immature male snow crab --", regions[i]), xlab="Year",
            ylab=expression(paste( "Number (X ", 10^6, " )")),
            cex.lab=cex.lab, cex.axis=cex.axis, cex.main = cex.main,
            panel = function(x, y, subscripts, ub, lb, ...) {
            larrows(x, lb[subscripts],
                    x, ub[subscripts],
                   angle = 90, code = 3, length=0.05)
            panel.xyplot(x, y, type="b", lty=1, lwd=2, pch=20, col="black", ...)
            panel.abline(v=2001.5, col="gray", lty=3, lwd=2, ...)
            panel.abline(h=0, col="black", lty=1, lwd=1, ...)
            }
          )
      print(pl)
    dev.off() 
    cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
 
    
    }
    return("Done")
  }  


