
  histograms.size.maturity.single.area = function( outdir, area, redo.data=F ) {
    # size frequency distributions of snow crab, broken down by maturity

      loc = file.path(project.datadirectory("snowcrab"), "R", "size.data")
       
      dir.create(path=outdir, recursive=T, showWarnings=F)
      dir.create(path=loc, recursive=T, showWarnings=F)
      outfilename = paste( c("mi", "mm", "fi", "fm"), "rdata", sep=".")
      outfile = file.path(loc, paste(outfilename))
      
      areas = area
      years = 1998:p$current.assessment.year
if(length(years)>15) years = (p$current.assessment.year-14):p$current.assessment.year
if(grepl('4x',tolower(areas))) years = 2004:p$current.assessment.year


      if (redo.data) {
      
        set = snowcrab.db( DS="set.clean")
        set$sid = paste(set$trip, set$set, sep="~")
     
        det = snowcrab.db( DS="det.initial")
        det$sid = paste(det$trip, det$set, sep="~")

        hvar="cw"
        bw = 2
        s0 = 10
        s1 = 150

        det=det[(det[,hvar]>=s0 & det[,hvar]<=s1),]
        breaks = seq(s0, s1, bw)

        m.imm = make.histograms(set, det[filter.class(det, "m.imm"),], hvar=hvar, breaks=breaks )
        m.mat = make.histograms(set, det[filter.class(det, "m.mat"),], hvar=hvar, breaks=breaks )
        f.imm = make.histograms(set, det[filter.class(det, "f.imm"),], hvar=hvar, breaks=breaks )
        f.mat = make.histograms(set, det[filter.class(det, "f.mat"),], hvar=hvar, breaks=breaks )

        save(m.imm, file=outfile[1], compress=T)
        save(m.mat, file=outfile[2], compress=T)
        save(f.imm, file=outfile[3], compress=T)
        save(f.mat, file=outfile[4], compress=T)
      
      }
      
      set = snowcrab.db( DS="set.clean")
      set$sid = paste(set$trip, set$set, sep="~")
      for (f in  outfile) load(f)


      # males
      fn = file.path(  outdir, areas )
      
      Cairo( file=fn, type="pdf", bg="white", units="in", width=8, height=10 )
      sexes = c(0,1)
      ncols = length(sexes)
      nrows = length(years)
      pl = layout( matrix( c(1:(ncols*nrows)), nrow=nrows, ncol=ncols, byrow=F ) )
      par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
      par(mar=c(0.4, 0, 0.4, 1.2))

      # ylim=c(0,400) # for 4X
      cols = c("gray40", "gray100" )
        set0 = set[filter.region.polygon(set, areas),]
      for (a in 1:(ncols)) {
    if(a==1){ ylim=c(0,750); xlim=c(0,140)}
    if(a==2){ ylim=c(0,1800); xlim=c(10,84)}
        for (y in 1:nrows) {
          set1 = set0[ which(set0$yr==years[y] ), ]
          sids = sort(unique(set1$sid))
          mim = f.imm; mm = f.mat
          if(sexes[a] == 0 ) {mim = m.imm; mm = m.mat}
          m.i = mim[which( rownames(mim)%in% sids ) ,]
          m.i.means = apply(X=m.i, MARGIN=2, FUN=mean, na.rm=T)
          m.m = mm[which( rownames(mm)%in% sids ) ,]
          m.m.means = apply(X=m.m, MARGIN=2, FUN=mean, na.rm=T)

          toplot = rbind(m.m.means, m.i.means)
          rn = as.numeric(colnames(toplot))
          toplot = toplot[, rn>=xlim[1] & rn<=xlim[2]]

          axes = T
          
          axisnames = F
          if (years[y]==years[nrows]) axisnames=T  # last row

          barplot(toplot, space=0, axisnames=axisnames, ylim=ylim, axes=axes, col=cols, xpd=F, lwd=3)
          
          if (sexes[a]==sexes[ncols]) {
            text( dim(toplot)[2]-4, ylim[2]*2/3, years[y], cex=1.2  )
          }
          
          if (sexes[a]==sexes[2] & years[y] %in% c(1998:2000) ) {
          } else {
            abline( v=41, lwd=3, lty="dashed" )
          }
          
          if (sexes[a]==sexes[2] & years[y]==years[1] ) {
            xl = c(xlim[2]*0.1, xlim[2]*0.1)
            yl = c(ylim[2]*0.8, ylim[2]*0.4 )
            points( x=xl, y=yl, pch=22, bg=c(cols[2], cols[1]), cex=2 )
            text( x=xl+xlim[2]*0.02, y=yl-ylim[2]*0.05, c("Immature", "Mature"), cex=1, pos=4)
          }

        }
     }

     mtext("Carapace width (mm)", side=1, outer=T, line=4, cex=1.2)
     mtext(expression(paste("No. / ", km^2)), side=2, outer=T, line=4, cex=1.2)
     mtext("Male", side=3, outer=T, line=1, at=0.33, cex=1.2)
     mtext("Female", side=3, outer=T, line=1, at=0.66, cex=1.2)
     
  dev.off()
  cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
   
  return("Done")
  }



