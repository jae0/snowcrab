
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  
  det = snowcrab.db("det.georef")

  hvar="cw" #!

  base = 2
  s0 = 10  # ~ ln(10)
  s1 = 160  # ~ ln(160)

  det=det[(det[,hvar]>=s0 & det[,hvar]<=s1),]

  s0 = log(s0, base=base)
  s1 = log(s1, base=base)
  det$cw = log(det$cw, base=base)

  bwh = 0.01 # bw for histogram
  bwd = 0.025 # bw for density

  m.tout = NULL
  f.tout = NULL

  years=c(1996:2006)
  areas = c("cfaall", "cfanorth", "cfa23", "cfa24", "cfa4x", "cfaslope")
  det0 = det

  x11()
  mm = which(det0$sex==1)
  mmh = jae.hist( det0$cw[mm], s0, s1, bwh, bwd, title="Male all",
    xlab=paste("CW (log base ", base, ")", sep=""))

  x11()
  for (y in years) {
  for (a in areas) {
  for (mat in c(1,2) ) {
    ij = filter.region.polygon(det0, a)
    if (length(ij)<30) next
    det = det0[ij,]
    ii = which(det$sex==1 & det$mat==mat & as.numeric(substring(det$trip,6,9))==y)
    if (length(ii)<30) next

    xxh = jae.hist( det$cw[ii], s0, s1, bwh, bwd, title="", xlab=paste("CW (log base ", base, ")", sep=""))
    print ("troughs")
    m.troughs = NULL
    m.troughs = locator()
    print(m.troughs)
    m.troughs$yr = y
    m.troughs$area = a
    m.tout = rbind(m.tout, as.data.frame(m.troughs))

  }}}
  save(m.tout, file="male.troughs.rdata")  ########### saved file!


  g = density(m.tout$x,bw=0.04)
  plot(g, col="red")
  male.res = locator()
  save(male.res, file="male.res.rdata") ########### saved file!

  male.res$x = c(3.232735,3.601103,4.025156,4.483475,4.946078,5.400114, 5.751349,7.194841)
  male.res$y = c(0.3135413, 0.6314096,0.7020470, 0.5557267, 0.6415007, 0.6667283, 0.5052714, 0.8433219)

 # h = base^male.res$x [ c(1:7)]  # <-----convert back to original units (cw; mm)
 # i = c(1:7)
 # the last point (index=8) is ambiguous ... upper bound is indeterminate

# based upon kernal density analysis and subsequent regression of
# size-class composition break-points in each of the areas and by time,
# optimal break points were identified (1996-2006) via the following regression:
#
# lm(formula = log(h) ~ i)  # <<<<------------- now using natural logarithm
# #lm(formula = log(cw;mm) ~ index(1:10)
# #  where index =1 was the first reliably identifiable left bound
#
#
# Residuals:
#         1         2         3         4         5         6         7
#  0.024282 -0.019298 -0.024282 -0.005514  0.016223  0.032023 -0.023434
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) 1.917564   0.022103   86.75 3.86e-09 ***
# i           0.298914   0.004942   60.48 2.34e-08 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
# Residual standard error: 0.02615 on 5 degrees of freedom
# Multiple R-Squared: 0.9986,     Adjusted R-squared: 0.9984
# F-statistic:  3658 on 1 and 5 DF,  p-value: 2.339e-08
#

# this function has been [;a;ced in "functions.trawl.r
#
 mb = function(i) exp(1.917564 + 0.298914*i)
    # mb = male bound in cw where i=1 when cw = ~ 9.2 mm,
    # instar 5 is bounded by mb(2) and mb(3)
    # prior to instar 5 moulting can be variable and so should not be
    # relied upon to heavily





# --------------------------------
# --- females


  det = snowcrab.db("det.georef")

  hvar="cw" #!

  base = 2
  s0 = 10  # ~ ln(10)
  s1 = 160  # ~ ln(160)

  det=det[(det[,hvar]>=s0 & det[,hvar]<=s1),]

  s0 = log(s0, base=base)
  s1 = log(s1, base=base)
  det$cw = log(det$cw, base=base)

  bwh = 0.01 # bw for histogram
  bwd = 0.025 # bw for density

  m.tout = NULL
  f.tout = NULL

  years=c(1996:2006)
  areas = c("cfaall", "cfanorth", "cfa23", "cfa24", "cfa4x", "cfaslope")
  det0 = det

  x11()
  mm = which(det0$sex==2)
  mmh = jae.hist( det0$cw[mm], s0, s1, bwh, bwd, title="Female all",
    xlab=paste("CW (log base ", base, ")", sep=""))

  x11()
  for (y in years) {
  for (a in areas) {
  for (mat in c(1,2) ) {
    ij = filter.region.polygon(det0, a)
    if (length(ij)<30) next
    det = det0[ij,]
    ii = which(det$sex==2 & det$mat==mat & as.numeric(substring(det$trip,6,9))==y)
    if (length(ii)<30) next

    xxh = jae.hist( det$cw[ii], s0, s1, bwh, bwd, title="", xlab=paste("CW (log base ", base, ")", sep=""))
    print ("troughs")
    f.troughs = NULL
    f.troughs = locator()
    print(f.troughs)
    f.troughs$yr = y
    f.troughs$area = a
    f.tout = rbind(f.tout, as.data.frame(f.troughs))

  }}}
  save(f.tout, file="female.troughs.rdata")  ########### saved file!


  g = density(f.tout$x,bw=0.04)
  plot(g, col="red")
  female.res = locator()
  save(female.res, file="female.res.rdata") ########### saved file!

  female.res$x = c( 3.642685, 3.825175, 4.058356, 4.524719, 5.011358, 5.447305, 5.893391, 6.359753)
  female.res$y = c(0.3515191, 0.4496409, 0.6090888, 1.0383715, 0.7440062, 1.2346151, 0.5354975, 0.4005800)

 # h = base^female.res$x [ c(1,3:7)]
# i = c(1:6)
 # -- the second point may be an additionsl moult within the earliest size classes (annual moults seem to begin at index=3)
 # the last point (index=8) is ambiguous ... upper bound is indeterminate
 # model = lm(log(h) ~ i)
# lm(formula = log(h) ~ i)
#
# Residuals:
#         1         2         3         4         5         6
#  0.011043 -0.015862 -0.007630  0.014657  0.001807 -0.004016
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) 2.198848   0.012012   183.1 5.34e-09 ***
# i           0.315026   0.003084   102.1 5.51e-08 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
# Residual standard error: 0.0129 on 4 degrees of freedom
# Multiple R-Squared: 0.9996,     Adjusted R-squared: 0.9995
# F-statistic: 1.043e+04 on 1 and 4 DF,  p-value: 5.511e-08
#

# this function has been [;a;ced in "functions.trawl.r
   fb = function(i) exp(2.198848 + 0.315026*i)
    # fb = female bound in cw where i=1 when cw = ~ 12.5 mm,
    # instar 5 is bounded by fb(1) and fb(2)
    # prior to instar 5 moulting can be variable and so should not be
    # relied upon to heavily






  m.tout$sex = "male"
  f.tout$sex = "female"

 m.tout$class = "trough"
 f.tout$class = "trough"


  sstruct = rbind(m.pout, m.tout, f.pout, f.tout)

  hist(det$cw[m1], breaks=breaks)
  hist(det$cw[m2], breaks=breaks)

  set$sid = paste(set$trip, set$set, sep="~")
  det$sid = paste(det$trip, det$set, sep="~")


  m.imm = make.histograms(set, det[filter.class(det, "m.imm"),], hvar=hvar, breaks=breaks )
  m.mat = make.histograms(set, det[filter.class(det, "m.mat"),], hvar=hvar, breaks=breaks )
  f.imm = make.histograms(set, det[filter.class(det, "f.imm"),], hvar=hvar, breaks=breaks )
  f.mat = make.histograms(set, det[filter.class(det, "f.mat"),], hvar=hvar, breaks=breaks )


  areas = "cfasouth"
  years = sort(unique(set$yr))


      # --------------------------
      # males

      ncols = length(areas)
      nrows = length(years)
      pl = layout( matrix( c(1:(ncols*nrows)), nrow=nrows, ncol=ncols, byrow=F ) )
      par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
      par(mar=c(0, 0, 0.4, 0))

      # ylim=c(0,400) # for 4X
      ylim=c(0,500)
      xlim=c(s0,s1)

      for (a in 1:(ncols)) {
        set0 = set[filter.region.polygon(set, areas[a]),]
        for (y in 1:nrows) {
          set1 = set0[ which(set0$yr==years[y]) , ]
          sids = sort(unique(set1$sid))

          m.i = m.imm[which( rownames(m.imm)%in% sids ) ,]
          m.i.means = apply(X=m.i, MARGIN=2, FUN=mean)
          m.m = m.mat[which( rownames(m.mat)%in% sids ) ,]
          m.m.means = apply(X=m.m, MARGIN=2, FUN=mean)

          toplot = rbind(m.m.means, m.i.means)
          rn = as.numeric(colnames(toplot))
          toplot = toplot[, rn>=xlim[1] & rn<=xlim[2]]

          axes = F
          if (areas[a]==areas[1] ) axes=T  # first col

          axisnames = F
          if (years[y]==years[nrows]) axisnames=T  # last row

          ltext = NULL
          if (areas[a]==areas[1] & years[y]==years[1] ) ltext = c("Mature", "Immature")

          barplot(toplot, space=0, axisnames=axisnames, ylim=ylim, axes=axes, legend=ltext, bty="n", xpd=F)

          if (areas[a]==areas[ncols]) text( dim(toplot)[2]-4, ylim[2]*2/3, years[y], cex=1.2 )
          }
        }

     mtext("Carapace width (mm)", side=1, outer=T, line=4, cex=1.2)
     mtext(expression(paste("No. / ", km^2)), side=2, outer=T, line=4, cex=1.2)

     mtext("Scotian Shelf", side=3, outer=T, line=1, at=0.15)
     mtext("Male", side=3, outer=T, line=4, cex=1.5)




#      Pr("png", dname="size.freq", fname="male", width=8, height=10)
#     Pr("png", dname="size.freq", fname="male.gb", width=8, height=10)



      # --------------------------
      # females
      ncols = length(areas)
      nrows = length(years)
      pl = layout( matrix( c(1:(ncols*nrows)), nrow=nrows, ncol=ncols, byrow=F ) )
      par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
      par(mar=c(0, 0, 0.4, 0))

      ylim=c(0,500)
      xlim=c(s0,s1)

      for (a in 1:(ncols)) {
        set0 = set[filter.region.polygon(set, areas[a]),]
        for (y in 1:nrows) {
          set1 = set0[ which(set0$yr==years[y]) , ]
          sids = sort(unique(set1$sid))

          f.i = f.imm[which( rownames(f.imm)%in% sids ) ,]
          f.i.means = apply(X=f.i, MARGIN=2, FUN=mean)
          f.m = f.mat[which( rownames(f.mat)%in% sids ) ,]
          f.m.means = apply(X=f.m, MARGIN=2, FUN=mean)

          toplot = rbind(f.m.means, f.i.means)
          rn = as.numeric(colnames(toplot))
          toplot = toplot[, rn>=xlim[1] & rn<=xlim[2]]

          axes = F
          if (areas[a]==areas[1] ) axes=T  # first col

          axisnames = F
          if (years[y]==years[nrows]) axisnames=T  # last row

          ltext = NULL
          if (areas[a]==areas[1] & years[y]==years[1] ) ltext = c("Mature", "Immature")

          barplot(toplot, space=0, axisnames=axisnames, ylim=ylim, axes=axes, legend=ltext, bty="n", xpd=F)

          if (areas[a]==areas[ncols]) text( dim(toplot)[2]-4, ylim[2]*2/3, years[y], cex=1.2 )
     }}

     mtext("Carapace width (mm)", side=1, outer=T, line=4, cex=1.2)
     mtext(expression(paste("No. / ", km^2)), side=2, outer=T, line=4, cex=1.2)

     mtext("Scotian Shelf", side=3, outer=T, line=1, at=0.15)
     mtext("Female", side=3, outer=T, line=4, cex=1.5)

#     Pr("png", dname="size.freq", fname="female", width=8, height=10)
   # Pr("png", dname="size.freq", fname="female.gb", width=8, height=10)





  o = r = NULL
  k = 1

  for (y in 1996:2006) {

    mi = which(det$sex==1 & det$mat==2 & as.numeric(substring(det$trip,6,9))==y)
    o = density(det$cw[mi], bw=0.02)
    o$y = o$y * k
    if (y==1996) {
      r = o
    }
    else {
      r$y = r$y + o$y
    }
  }

  plot(r)

