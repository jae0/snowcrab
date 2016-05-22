
  figure.growth.instar = function( outdir) {
    
    instars = c(2:12)  # 3 is instar 6
    y = mb(instars)
    ym = (y[2:12] + y[1:11]) / 2

    ym=ym[ym<=170]
    ym = ym[is.finite(ym)]

    growth = data.frame( instar=c(5:13), cw=ym)
    growth$mass = growth$cw^ 3.098 *  0.0002665
    growth$mass[ which(growth$instar==13)] = growth$mass[ which(growth$instar==12)] + 450  # the last increment does not seem allometric .. computed from empirical data of CC13-CC12
    growth$mass = growth$mass / 1000

    fn = file.path( outdir, "male.growth.cw" )
    Cairo( file=fn, type="pdf", bg="white", units="in", width=8, height=6 )
    plot( growth$instar, growth$cw, type="b", ylab="Carapace width (mm)", xlab="Instar")
    dev.off()
    cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
   
    fn = file.path( outdir, "male.growth.mass" )
    Cairo( file=fn, type="pdf", bg="white", units="in", width=8, height=6 )
    plot( growth$instar, growth$mass, type="b", ylab="Body mass (kg)", xlab="Instar")
    dev.off()
    cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
 
  }


