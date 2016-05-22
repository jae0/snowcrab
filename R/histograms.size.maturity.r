
  histograms.size.maturity = function( outdir, redo.data=F ) {
    # size frequency distributions of snow crab, broken down by maturity

      loc = file.path(project.datadirectory("snowcrab"), "R", "size.data")
       
      dir.create(path=outdir, recursive=T, showWarnings=F)
      dir.create(path=loc, recursive=T, showWarnings=F)
      outfilename = paste( c("mi", "mm", "fi", "fm"), "rdata", sep=".")
      outfile = file.path(loc, paste(outfilename))
       
      # areas = c("cfanorth.not.glace.bay", "cfa22outer", "cfasouth" )
      # areas = c("cfaall",  "cfanorth", "cfasouth", "cfa20", "cfa21", "cfa22", "cfa23", "cfa24", "cfa4x", "cfa23slope", "cfa24slope", "cfaslope"  )
      areas = c("cfanorth", "cfasouth", "cfa4x")
      year = 1998:p$current.assessment.year
    if(length(year)>15) year = (p$current.assessment.year-14):p$current.assessment.year



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
      fn = file.path(  outdir, "male" )
      
      Cairo( file=fn, type="pdf", bg="white", units="in", width=8, height=10 )

      ncols = length(areas)
      nrows = length(year)
      pl = layout( matrix( c(1:(ncols*nrows)), nrow=nrows, ncol=ncols, byrow=F ) )
      par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
      par(mar=c(0, 0, 0.4, 0))

      # ylim=c(0,400) # for 4X
      ylim=c(0,750)
      xlim=c(0,140)
      cols = c("gray40", "gray100" )

      for (a in 1:(ncols)) {
        set0 = set[filter.region.polygon(set, areas[a]),]
        for (y in 1:nrows) {
          set1 = set0[ which(set0$yr==year[y] ), ]
          sids = sort(unique(set1$sid))

          m.i = m.imm[which( rownames(m.imm)%in% sids ) ,]
          m.i.means = apply(X=m.i, MARGIN=2, FUN=mean, na.rm=T)
          m.m = m.mat[which( rownames(m.mat)%in% sids ) ,]
          m.m.means = apply(X=m.m, MARGIN=2, FUN=mean, na.rm=T)

          toplot = rbind(m.m.means, m.i.means)
          rn = as.numeric(colnames(toplot))
          toplot = toplot[, rn>=xlim[1] & rn<=xlim[2]]

          axes = F
          if (areas[a]==areas[1] ) axes=T  # first col

          axisnames = F
          if (year[y]==year[nrows]) axisnames=T  # last row

          barplot(toplot, space=0, axisnames=axisnames, ylim=ylim, axes=axes, col=cols, xpd=F, lwd=3)
          
          if (areas[a]==areas[ncols]) {
            text( dim(toplot)[2]-4, ylim[2]*2/3, year[y], cex=1.2  )
          }
          
          if (areas[a]==areas[3] & year[y] %in% c(1998:2000) ) {
          } else {
            abline( v=41, lwd=3, lty="dashed" )
          }
          
          if (areas[a]==areas[3] & year[y]==year[2] ) {
            xl = c(xlim[2]*0.1, xlim[2]*0.1)
            yl = c(ylim[2]*0.8, ylim[2]*0.4 )
            points( x=xl, y=yl, pch=22, bg=c(cols[2], cols[1]), cex=2 )
            text( x=xl+xlim[2]*0.02, y=yl-ylim[2]*0.05, c("Immature", "Mature"), cex=1, pos=4)
          }

        }
     }

     mtext("Carapace width (mm)", side=1, outer=T, line=4, cex=1.2)
     mtext(expression(paste("No. / ", km^2)), side=2, outer=T, line=4, cex=1.2)
     mtext("N-ENS", side=3, outer=T, line=1, at=0.15, cex=1.2)
     mtext("S-ENS", side=3, outer=T, line=1, at=0.5, cex=1.2)
     mtext("4X", side=3, outer=T, line=1, at=0.85, cex=1.2)
     mtext("MALE", side=3, outer=T, line=4, cex=1.4)

  dev.off()
  cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
   
      # females
  fn = file.path(  outdir, "female" )
  Cairo( file=fn, type="pdf", bg="white", units="in", width=8, height=10 )
    
      ncols = length(areas)
      nrows = length(year)
      pl = layout( matrix( c(1:(ncols*nrows)), nrow=nrows, ncol=ncols, byrow=F ) )
      par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
      par(mar=c(0, 0, 0.4, 0))

      # ylim=c(0,400)
      ylim=c(0,1800)
      xlim=c(10,84)

      cols = c("gray40", "gray100" )

      for (a in 1:(ncols)) {
        set0 = set[filter.region.polygon(set, areas[a]),]
        for (y in 1:nrows) {
          set1 = set0[ which(set0$yr==year[y]) , ]
          sids = sort(unique(set1$sid))

          f.i = f.imm[which( rownames(f.imm)%in% sids ) ,]
          f.i.means = apply(X=f.i, MARGIN=2, FUN=mean, na.rm=T)
          f.m = f.mat[which( rownames(f.mat)%in% sids ) ,]
          f.m.means = apply(X=f.m, MARGIN=2, FUN=mean, na.rm=T)

          toplot = rbind(f.m.means, f.i.means)
          rn = as.numeric(colnames(toplot))
          toplot = toplot[, rn>=xlim[1] & rn<=xlim[2]]

          axes = F
          if (areas[a]==areas[1] ) axes=T  # first col

          axisnames = F
          if (year[y]==year[nrows]) axisnames=T  # last row

          barplot(toplot, space=0, axisnames=axisnames, ylim=ylim, axes=axes, xpd=F, lwd=3 )

          if (areas[a]==areas[ncols]) text( dim(toplot)[2]-4, ylim[2]*2/3, year[y], cex=1.2 )

          if (areas[a]==areas[3] & year[y]==year[2] ) {
            xl = c(xlim[2]*0.1, xlim[2]*0.1)
            yl = c(ylim[2]*0.8, ylim[2]*0.4 )
            points( x=xl, y=yl, pch=22, bg=c(cols[2], cols[1]), cex=1 )
            text( x=xl+xlim[2]*0.02, y=yl-ylim[2]*0.05, c("Immature", "Mature"), cex=1, pos=4)
          }

     }}

      mtext("Carapace width (mm)", side=1, outer=T, line=4, cex=1.2)
      mtext(expression(paste("No. / ", km^2)), side=2, outer=T, line=4, cex=1.2)
      mtext("N-ENS", side=3, outer=T, line=1, at=0.15, cex=1.2)
      mtext("S-ENS", side=3, outer=T, line=1, at=0.5, cex=1.2)
      mtext("4X",   side=3, outer=T, line=1, at=0.85, cex=1.2)
      mtext("FEMALE", side=3, outer=T, line=4, cex=1.2)
      
    dev.off()
 cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
 
  return("Done")
  }



