

	
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
 

  # white hake
  wh2 = snowcrab.db( DS="cat.georeferenced" )
  wh2$id = paste(wh2$trip, wh2$set, wh2$yr )
  
  wh = wh2[ -which(duplicated(wh2$id)) , c( "id", "lon", "lat", "yr" )]
  
  wh2 = wh2[ which (wh2$spec == 12 ) , ]
  wh = merge( wh, wh2[, c("id", "totno", "totmass") ], by="id", all.x=T, all.y=F )
  wh = wh[ which(wh$yr>= 2004 ) ,]

  write.csv( wh, file="~/whitehake.csv")



 # Fred Kennedy's request
 # Carapace condition breakdown by year
  set = snowcrab.db("set.with.det")
  x = set[ filter.region(set$cfa, "cfa23"), ]

  vars = c("totmass.male.com.CC1", "totmass.male.com.CC2", "totmass.male.com.CC3", "totmass.male.com.CC4", "totmass.male.com.CC5")
  out =NULL
  for (y in unique(x$yr)) {
    a = which(x$yr == y)
    for (v in vars) {
      q = x[a,v]+1
      mcc = exp(mean(log(q), na.rm=T)) - 1
      lcc = exp(mcc - (sd  (log(q), na.rm=T))) - 1
      ucc = exp(mcc + (sd  (log(q), na.rm=T))) - 1
      ncc = length(q)
      cc  = cbind(y, v, mcc, lcc, ucc, ncc)
      out = rbind(out, cc)
    }}

  write.table(out, file="tmp.csv", quote=F, col.names=T)


  # ----------------------------------------------------------
  # Roger Pettipas .. bottom temperatures of Scotian shelf
  set = snowcrab.db( "set.merge.det" )
  varstoextract = c("t0", "lon", "lat", "z", "zsd", "t", "tsd")
  t = set[ which(set$yr==2012), varstoextract]

  write.table(t, file="petipas.csv", col.names = T, row.names=F, quote=F, sep=" ; ")


  # ----------------------------------------------------------
  # Peter Koeller -- N immature male > 56 mm CW in areas 23ab and 24ab
  #  == R2 + R3 + R4
     
  loadfunctions( "snowcrab", functionname="initialise.local.environment.r")
  
  set = snowcrab.db("set.complete")
  # overrides:
  p$regions.to.model = "cfa.23ab.24ab"
  p$vars.to.model = "pre.recruit.no"
  p$years.to.model = p$years.to.model[ which(p$years.to.model>1997) ]
  p$ofname = "kriged.results.koeller.rdata"

  pcoords = c("plon", "plat")

  # new method: directly computed averages of core areas
  i = filter.region.polygon(x=set[, pcoords], region=p$regions.to.model, planar=T, proj.type=p$internal.projection )
  xs = set[ i, ]
  xs$dummy = 1

  res = tapply( xs$pre.recruit.no, INDEX=xs$yr, FUN=mean, na.rm=T ) 
  res2 = tapply( xs$pre.recruit.no, INDEX=xs$yr, FUN=sd, na.rm=T ) 
  res3 = tapply( xs$dummy, INDEX=xs$yr, FUN=sum, na.rm=T ) 
   
  res = as.data.frame( cbind(res, res2, res3 ) )
  names(res) = c("mean", "sd", "n")
  res$se = res$sd/ sqrt(res$n+1)

  out = res[, c("mean", "se")]
  out$year = as.numeric(rownames(out))
  rownames(out)=NULL
  write.table ( out, file="~/tmp/snowcrab.recruitment.index.csv", sep="\t" )
      
  plot(  res[, c("mean")], ylim=c(0, max(res$mean, na.rm=T)*1.1))


  # ---------------------------
  # Glace Bay Hole ... briefing note

  # for Glace Bay comparisons: override with ..
  
  # reulsts in: ~/ecomod/snowcrab/issues/briefing.note.glace.bay/2005
  
  loadfunctions( "snowcrab", functionname="initialise.local.environment.r")
  
  set = snowcrab.db("set")

  # overrides:
  p$regions.to.model = c( "cfanorth", "cfasouth", "cfa22outer", "cfa23a", "cfanorth.not.glace.bay" )
  p$vars.to.model = c( "totmass.male.com", "R0.no", "totno.female.mat", "totno.female.berried", "totno.all")
  p$years.to.model = p$years.to.model[which(p$years.to.model>1997)]
  p$ofname = "kriged.results.glacebay2.rdata"
  p$krige.ordinary = F
  p$krige.block = T
  p$clusters=rep("tethys",4)  # using all 8 will make it run out of memory
 
  K = krige.block ( p ) # krige ordinary should already have been done


  for (v in p$vars.to.model) {
  for (r in p$regions.to.model) {
    print (paste(v,r))
    XX=K[ which(K$vars==v & K$region==r) ,]
    xx = XX$yr
    yy = XX$total
    lbound = XX$lbound
    ubound = XX$ubound
    surfacearea = XX$surfacearea

    ts.plotandsave(xx, yy, lbound, ubound, surfacearea, action="save", title="Biomass (x 1000 t)", filename=paste(v,r,sep="."), outdir=file.path(project.datadirectory("snowcrab"), "R", "Requests") )

  }
  }

  # ------------------------
  # size freq distributions

	loadfunctions( "snowcrab", functionname="initialise.local.environment.r")
     
  loc = file.path(p$annual.results, "size.data")
 
  dir.create(path=loc, recursive=T, showWarnings=F)
  outfilename = paste( c("mi", "mm", "fi", "fm"), "rdata", sep=".")
  outfile = file.path(loc, paste(outfilename))

  if (!redo.data) {
     
     for (f in  outfile) load(f)

  } else if (redo.data) {

    det = snowcrab.db("det")
    set = snowcrab.db("set.netmind.sanity.checked")
    v = c( "male", "female")
    varstokeep = c("lon", "lat", "trip", "set", "yr", "cw", "sa.x", "sex", "shell", "durometer", "mat", "stage.fmat", "eggPr")

      # Glace Bay extraction
      areas = c("cfanorth.not.glace.bay", "cfa22outer", "cfasouth" )
      years = 1998:p$current.assessment.year

      set$sid = paste(set$trip, set$set, sep="~")
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
  # --------------------------
  # males

  ncols = length(areas)
  nrows = length(years)
  pl = layout( matrix( c(1:(ncols*nrows)), nrow=nrows, ncol=ncols, byrow=F ) )
  par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
  par(mar=c(0, 0, 0.4, 0))

  # ylim=c(0,400) # for 4X
  ylim=c(0,1200)
  xlim=c(0,140)
   

  cols = c("gray40", "gray100" )

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

      barplot(toplot, space=0, axisnames=axisnames, ylim=ylim, axes=axes, col=cols, xpd=F)

      if (areas[a]==areas[ncols]) {
        text( dim(toplot)[2]-4, ylim[2]*2/3, years[y], cex=1.2 )
      }

      if (areas[a]==areas[3] & years[y]==years[2] ) {
        xl = c(xlim[2]*0.1, xlim[2]*0.1)
        yl = c(ylim[2]*0.8, ylim[2]*0.4 )
        points( x=xl, y=yl, pch=22, bg=c(cols[2], cols[1]), cex=2 )
        text( x=xl+xlim[2]*0.02, y=yl-ylim[2]*0.05, c("Immature", "Mature"), cex=1, pos=4)
      }

    }
 }

 mtext("Carapace width (mm)", side=1, outer=T, line=4, cex=1.2)
 mtext(expression(paste("No. / ", km^2)), side=2, outer=T, line=4, cex=1.2)
 mtext("North", side=3, outer=T, line=1, at=0.15, cex=1.2)
 mtext("GBH",   side=3, outer=T, line=1, at=0.5, cex=1.2)
 mtext("South", side=3, outer=T, line=1, at=0.85, cex=1.2)
 mtext("MALE", side=3, outer=T, line=4, cex=1.4)

 Pr("png", dname="size.freq", fname="male.gb", width=8, height=10)



  # --------------------------
  # females
  
  ncols = length(areas)
  nrows = length(years)
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

      barplot(toplot, space=0, axisnames=axisnames, ylim=ylim, axes=axes, xpd=F)

      if (areas[a]==areas[ncols]) text( dim(toplot)[2]-4, ylim[2]*2/3, years[y], cex=1.2 )

      if (areas[a]==areas[3] & years[y]==years[2] ) {
        xl = c(xlim[2]*0.1, xlim[2]*0.1)
        yl = c(ylim[2]*0.8, ylim[2]*0.4 )
        points( x=xl, y=yl, pch=22, bg=c(cols[2], cols[1]), cex=1 )
        text( x=xl+xlim[2]*0.02, y=yl-ylim[2]*0.05, c("Immature", "Mature"), cex=1, pos=4)
      }

 }}

 mtext("Carapace width (mm)", side=1, outer=T, line=4, cex=1.2)
 mtext(expression(paste("No. / ", km^2)), side=2, outer=T, line=4, cex=1.2)
 mtext("North", side=3, outer=T, line=1, at=0.15, cex=1.2)
 mtext("GBH",   side=3, outer=T, line=1, at=0.5, cex=1.2)
 mtext("South", side=3, outer=T, line=1, at=0.85, cex=1.2)
 mtext("FEMALE", side=3, outer=T, line=4, cex=1.2)

 Pr("png", dname="size.freq", fname="female.gb", width=8, height=10)



  ### -- GBH continued ... movement
  # source the file "mark.recapture.r"

 k = move$dx / as.numeric(move$dt) * 365
 m = mean(log(k))
 sd = sqrt(var(log(k), na.rm=T))
 ub = exp( m + 2*sd )
 lb = exp( m - 2*sd )
 hist( log (k), "FD")

 x = move

  # starting locations from within glace bay ( "cfa22outer" )
  x$lon = x$lon0
  x$lat = x$lat0
  i = filter.region.polygon( x, region="cfa22outer" )

  x = move
  # recaps in glace bay ( "cfa22outer" )
  x$lon = x$lon1
  x$lat = x$lat1
  i = filter.region.polygon( x, region="cfa22outer" )


# ---------------------
# data request: winter skate in the commercial observer database
# data extraction from orable db directly

observer.data.request.oracle = function () {

    res = dbGetQuery(con, paste(
  "SELECT trip.trip_id, trip.trip, trip.board_date, st.set_no, ca.speccd_id, ca.est_num_caught, ca.est_kept_wt, ca.est_discard_wt, fish.fish_no , fish.fish_length, fish.fish_weight",
  "FROM istrips trip, isgears gr, isfishsets st, iscatches ca, isfish fish",
  "WHERE trip.tripcd_id = 2509",
  "AND trip.trip_id = gr.trip_Id",
  "AND (trip.trip_id = st.trip_Id AND gr.gear_id = st.gear_id)",
  "AND st.fishset_id = ca.fishset_id(+)",
  "AND ca.speccd_id(+) != 2526",
  "AND ca.catch_id = fish.catch_id(+)",
  "GROUP BY trip.trip_id, trip.trip, trip.board_date, st.set_no, ca.speccd_id, ca.est_num_caught, ca.est_kept_wt, ca.est_discard_wt, fish.fish_no , fish.fish_length, fish.fish_weight",
  "order by board_date, set_no, fish_no;"
          ) )
   f=res[ which(res$SPECCD_ID==204) ,]
   save(f, file="winter.skate.rdata", compress=T)
}


# ----------------------------------
# Hebert: mean size of mature females


		loadfunctions( "snowcrab", functionname="initialise.local.environment.r")

    det = snowcrab.db("det.georef")  # this contains year
    det = det[ which(det$mat ==1 & det$sex==2) ,]

    means = as.data.frame.table( tapply( det$cw, det$yr, mean, na.rm=T)  )
    names(means) = c("yr", "cw.mean")
    means$yr = as.numeric( as.character( means$yr) )

    sd = as.data.frame.table( tapply( det$cw, det$yr, sd, na.rm=T))
    names(sd) = c("yr", "cw.sd")
    sd$yr = as.numeric( as.character( sd$yr) )

    n = as.data.frame.table( tapply( det$cw, det$yr, function(x) {length( which(is.finite(x)))}) )
    names(n) = c("yr", "n")
    n$yr = as.numeric( as.character( n$yr) )

    res=NULL
    res = merge(means, sd, by="yr" ) 
    res = merge(res, n, by="yr" )
    res = res[res$yr > 1996 , ]
    
    plot( res$yr, res$cw.mean)



# ----------------------------------
# Gordon MacDonald: CFA 23/24 partitionning

	loadfunctions( "snowcrab", functionname="initialise.local.environment.r")

  set = snowcrab.db("set")

  # overrides:
  p$regions.to.model = c("cfa23", "cfa24")
  p$vars.to.model = "R0.mass"
  p$years.to.model = p$years.to.model[ which(p$years.to.model>2000) ]
  p$ofname = "kriged.results.CFA23.24.rdata"

  p$krige.ordinary = F # should already be done by this point
  p$krige.block = T  # if ordinary kriging already done ... to est CI

  p$do.parallel = F

  # p$clusters=rep("io",7)  # using all 8 will make it run out of memory
  
  if (p$krige.ordinary) krige.ordinary( p )  # ~ 6 hr on "io X 8"
  if (p$krige.block)    K = krige.block ( p ) # to est CI and get data summaries ~ 2 days

  # cmd("sync.jae io") 

  K = krige.block(p, source="file") # load 

  v = p$vars.to.model
  r = p$regions.to.model

  print( K[, c("vars", "region", "yr", "total", "lbound", "ubound")])
 
 
  cex.main = 1.4
  cex.lab = 1.3
  cex.axis = 1.3

  figure.height=8
  figure.width=6

  xlim=range(K$yr, na.rm=T)
  xlim[1]=xlim[1]-0.5
  xlim[2]=xlim[2]+0.5


  outdir = "ts.kriged.results"

  zero.value = K[1,]
  zero.value$yr = p$current.assessment.year
  zero.value$total = 0
  zero.value$lbound = 0
  zero.value$ubound = 0

  tt = "Fishable biomass"
  td = K
  
  yy = expression(paste( "Biomass (X ", 10^3, " t)"))
  convert = 1 # no conversion
  eps = 1e-6

  varstocheck = c("total", "ubound", "lbound")
  
  for (vs in varstocheck) {
        td[,vs] = td[,vs] / convert
        kk = which(td[,vs] <= eps)
        if (length(kk)>0) td[kk,vs] = 0
  }

  td = td[is.finite(td$total) ,]
  td = td[order(td$region,td$yr),]

      # last check to make sure data exists for all regions for the plots
      for (rr in 1:length(p$regions.to.model)) {
        nr=NULL
        nr = length( which(td$region==p$regions.to.model[rr] ))
        if (nr ==0) {
          mm = zero.value
          mm$region = p$regions.to.model[rr]
          td = rbind( td, mm)
        }
      }


      setup.lattice.options()
      y = xyplot( total~yr|region, data=td, ub=td$ubound, lb=td$lbound,
        layout=c(1,length(p$regions.to.model)), xlim=xlim, scales = list(y = "free"),
            main=tt, xlab="Survey year", ylab=yy, cex.lab=cex.lab, cex.axis=cex.axis, cex.main = cex.main,
            panel = function(x, y, subscripts, ub, lb, ...) {
            larrows(x, lb[subscripts],
                    x, ub[subscripts],
                   angle = 90, code = 3, length=0.05)
            panel.xyplot(x, y, type="b", lty=1, lwd=2, pch=20, col="black", ...)
            panel.abline(v=2001.5, col="gray", lty=1, lwd=2, ...)
            panel.abline(h=0, col="black", lty=1, lwd=1, ...)
            }
          )
     print(y)

     Pr(dev="png", dname=outdir, fname=v, width=figure.width, height=figure.height)
     options(default.options)

     l23 = td[ which(td$region =="cfa23"), ]
     l24 = td[ which(td$region =="cfa24"), ]
     l = merge(l23, l24, by="yr", suffixes=c(".23", ".24" ) )
      plot(l$yr, l$total.23/(l$total.24+l$total.23))



# ----------------------------------
# MacMullin: GBH vs NENS partitionning
 
  # results location: ~/ecomod/snowcrab/issues/briefing.note.glace.bay/2006
  
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r")

  set = snowcrab.db("set")

  # overrides:  
  p$regions.to.model = c( "cfanorth", "cfasouth", "cfa22outer", "cfa23a", "cfanorth.not.glace.bay" )
  p$vars.to.model = "R0.mass"
  p$years.to.model = p$years.to.model[ which(p$years.to.model>2000) ]
  p$ofname = "kriged.results.GBH.rdata"

  p$krige.ordinary = F # should already be done by this point
  p$krige.block = T  # if ordinary kriging already done ... to est CI

  p$do.parallel = F

  # p$clusters=rep("io",7)  # using all 8 will make it run out of memory
  
  if (p$krige.ordinary) krige.ordinary( p )  # ~ 6 hr on "io X 8"
  if (p$krige.block)    K = krige.block ( p ) # to est CI and get data summaries ~ 2 days

  # cmd("sync.jae io") 

  K = krige.block(p, source="file") # load 

  v = p$vars.to.model
  r = p$regions.to.model

  print( K[, c("vars", "region", "yr", "total", "lbound", "ubound")])
 
 
  cex.main = 1.4
  cex.lab = 1.3
  cex.axis = 1.3

  figure.height=8
  figure.width=6

  xlim=range(K$yr, na.rm=T)
  xlim[1]=xlim[1]-0.5
  xlim[2]=xlim[2]+0.5


  outdir = "ts.kriged.results"

  zero.value = K[1,]
  zero.value$yr = p$current.assessment.year
  zero.value$total = 0
  zero.value$lbound = 0
  zero.value$ubound = 0

  tt = "Fishable biomass"
  td = K
  
  yy = expression(paste( "Biomass (X ", 10^3, " t)"))
  convert = 1 # no conversion
  eps = 1e-6

  varstocheck = c("total", "ubound", "lbound")
  
  for (vs in varstocheck) {
        td[,vs] = td[,vs] / convert
        kk = which(td[,vs] <= eps)
        if (length(kk)>0) td[kk,vs] = 0
  }

  td = td[is.finite(td$total) ,]
  td = td[order(td$region,td$yr),]

      # last check to make sure data exists for all regions for the plots
      for (rr in 1:length(p$regions.to.model)) {
        nr=NULL
        nr = length( which(td$region==p$regions.to.model[rr] ))
        if (nr ==0) {
          mm = zero.value
          mm$region = p$regions.to.model[rr]
          td = rbind( td, mm)
        }
      }


      setup.lattice.options()
      y = xyplot( total~yr|region, data=td, ub=td$ubound, lb=td$lbound,
        layout=c(1,length(p$regions.to.model)), xlim=xlim, scales = list(y = "free"),
            main=tt, xlab="Survey year", ylab=yy, cex.lab=cex.lab, cex.axis=cex.axis, cex.main = cex.main,
            panel = function(x, y, subscripts, ub, lb, ...) {
            larrows(x, lb[subscripts],
                    x, ub[subscripts],
                   angle = 90, code = 3, length=0.05)
            panel.xyplot(x, y, type="b", lty=1, lwd=2, pch=20, col="black", ...)
            panel.abline(v=2001.5, col="gray", lty=1, lwd=2, ...)
            panel.abline(h=0, col="black", lty=1, lwd=1, ...)
            }
          )
     print(y)

     Pr(dev="png", dname=outdir, fname=v, width=figure.width, height=figure.height)
     options(default.options)



# Sherrylynn Rowe: sea cucumber and whelks in surveys:

	loadfunctions( "snowcrab", functionname="initialise.local.environment.r")

	load(file.path( project.datadirectory("snowcrab"), "R", "cat_georef.rdata"))
   
whelks = c(4210, 4211, 4212, 4227, 4228)  # Family Buccinidae 
cucumbers = c(6600, 6600, 6611, 6700 , 6705, 6710 , 6711 , 6712, 6713, 6714, 6715, 6716, 6717, 6718, 6719, 6720 ) # class Holothuroidea
to.extract = c("lon", "lat", "totno", "yr", "spec") 

W = cat[ cat$spec %in% whelks, to.extract ]
C = cat[ cat$spec %in% cucumbers, to.extract ]

W = W[ which(W$yr >= 2005) ,]
C = C[ which(C$yr >= 2005) ,]

W$spec = "Family.Buccinidae"
C$spec = "Class.Holothuroidea"

out = rbind(W,C)

write.table( out, file="temp.csv", sep =";")



# ---------------- Zwanenburg / MPAs etc
 

loadfunctions( "snowcrab", functionname="initialise.local.environment.r")
     
load(file.path( project.datadirectory("snowcrab"), "R", "cat.georef.rdata"))
cat$totno = cat$totno * cat$sa
to.extract = c("lon", "lat", "yr", "totno", "spec") 
yrs = which(cat$yr > 2003)
out = cat[ yrs, to.extract ]
write.table( out, file="temp.csv", sep =";")




# 4X trawl locations
#
loadfunctions( "snowcrab", functionname="initialise.local.environment.r")
     
  set = snowcrab.db("setInitial")
  ii = filter.region.polygon(set, region="cfa4x")
  set = set[ ii, ]
  set = set[ , c("yr", "lon", "lat") ]
  (set)



# wolf fish extraction for Jim Simon
loadfunctions( "snowcrab", functionname="initialise.local.environment.r")
  
Y = snowcrab.db( DS="cat.georeferenced" )
spec = taxonomy.recode( from="taxa", to="spec", tolookup="wolffish")
taxonomy.recode( from="spec", to="taxa", tolookup=spec )
                 
spec = c(50, 51, 52)

Y = Y[ which( Y$spec %in% spec & Y$yr %in% c(2009:2011) ) , ]
Y$totmass = NULL
Y$plon = NULL
Y$plat = NULL
Y$sa = NULL

write.csv( Y, file="~/tmp/jim.csv" )



---- 
# skate purse extraction for Jim Simon
loadfunctions( "snowcrab", functionname="initialise.local.environment.r")

yrs = c(2005:2012)
Z = snowcrab.db( DS="set.complete" )
Z = Z[which(Z$yr %in% yrs ) , ]
Z = Z[, c("trip", "set", "yr", "lon", "lat")]

Y = snowcrab.db( DS="cat.georeferenced" )
spec = c(1224)

Y = Y[ which( Y$spec %in% spec ) , ]

Y = Y[ which( Y$spec %in% spec & Y$yr %in% yrs ) , ]
Y = Y[, c("trip", "set", "totno")]

Y = merge( Z, Y, by=c("trip", "set"), all.x=T, all.Y=F )
Y$trip = NULL
Y$set = NULL
Y$totno[ which(!is.finite( Y$totno))] = 0


write.csv( Y, file="~/tmp/jim.csv" )


# ------
# NENS data request WRT EMERA power line

loadfunctions( "snowcrab", functionname="initialise.local.environment.r")
set = snowcrab.db("set")
set = set[,c( "yr", "lon", "lat", "t", "sa", "R0.mass", "totmass.female.mat", "totmass.male.imm", "totmass.female.imm" )]
nens = filter.region.polygon( set, region="cfanorth" )
innergutter = intersect( nens, which(set$lon < -59.8 ) )
nensdata = set[innergutter,]



