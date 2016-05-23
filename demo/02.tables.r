
  # Tables based upon data created by "1.snowcrab.r"
	
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  library("xtable")
  library("R2HTML")

  odb0 = observer.db("odb")
  regions = c("cfanorth", "cfasouth", "cfa4x")
  nregions = length(regions)

  #------------------------------------------------
  #Fisheries statistics per region
  tabledir = file.path(project.datadirectory("snowcrab"), "data", "fisheries")
  outtabledir= file.path(project.datadirectory("snowcrab"), "assessments", "2015", "tables", "logbook")
  setwd(tabledir)
  
  NFS <- xtable(read.csv("NENS_FisherySummary.csv"))
  SFS <- xtable(read.csv("SENS_FisherySummary.csv"))
  Fx <- xtable(read.csv("4x_FisherySummary.csv"))

  setwd(outtabledir)
  print.xtable(NFS, type="latex", file="NENS_FisherySummary.tex")
  print.xtable(NFS, type="html", file="NENS_FisherySummary.html")
  
  print.xtable(SFS, type="latex", file="SENS_FisherySummary.tex")
  print.xtable(SFS, type="html", file="SENS_FisherySummary.html")

  print.xtable(Fx, type="latex", file="4x_FisherySummary.tex")
  print.xtable(Fx, type="html", file="4x_FisherySummary.html")
  
  #regions = c("cfaall")
   #regions = c("cfanorth", "cfasouth", "cfa4x")
   regions = c("cfanorth", "cfa23", "cfa24", "cfa4x")
    l = NULL
    for (r in regions) {
      res = get.fishery.stats.by.region( Reg=r) #need to add the TACs per year and number of licences
      #round the landings to ton
      #round CPUE to no decimal places
      #round the effort to per x1000 trap hauls
      print(r)
      print(res)
      }
# ----------------------------------------
#  Carapace condition from observed data  < 95mm CW
   
    outtabledir= file.path(project.datadirectory("snowcrab"), "assessments", "2015", "tables", "observer")

    odb = odb0
    odb = odb[ which( odb$cw < 95 & odb$prodcd_id=="0" ) ,] 
    regions = c("cfanorth", "cfasouth", "cfa4x")
    nregions = length(regions) 
    years = sort( unique( odb$fishyr ) )

    res = NULL
    for (r in p$regions) {
    for (y in years) {
      out = proportion.cc (odb, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

    cnames = c("region", "fishyr", c(1:5), "ntot")
    colnames(res) = cnames
    print(res)
    res = as.data.frame(res)
    res[is.na(res)] <- NA
    ct <- c("CC1", "CC2", "CC3", "CC4", "CC5", "Total")

    setwd(outtabledir)

    Rn= res[res$region=="cfanorth", 3:8]
    print(Rn)
    rownames(Rn) = years
    colnames(Rn) = ct
    print.xtable(Rn, type="latex", file="table.CC.small.north.obs.tex")
    HTML(Rn, file="table.CC.Small.north.obs.html")

    Rs= res[res$region=="cfasouth", 3:8]
    rownames(Rs) = years
    colnames(Rs) = ct
    print.xtable(Rs, type="latex", file="table.CC.small.south.obs.tex")
    HTML(Rs, file="table.CC.small.south.obs.html")

    Rx= res[res$region=="cfa4x", 3:8]
    rownames(Rx) = years
    colnames(Rx) = ct
    print.xtable(Rs, type="latex", file="table.CC.small.4x.obs.tex")
    HTML(Rx, file="table.CC.small.4x.obs.html")
# ----------------------------------------
#  Carapace condition from observed data >=95mm CW
    odb = odb0
    odb = odb[ which( odb$cw >= 95 & odb$cw < 170 & odb$prodcd_id=="0" ) ,]  # commerical sized crab only
    years = sort( unique( odb$fishyr ) )

  # get proportion by cc
    regions = c("cfanorth", "cfasouth", "cfa4x")
    years = sort( unique( odb$fishyr ) )

    res = NULL
    for (r in regions) {
    for (y in years) {
      out = proportion.cc (odb, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

    cnames = c("region", "fishyr", c(1:5), "ntot")
    colnames(res) = cnames
    print(res)
    res = as.data.frame(res)
    res[is.na(res)] <- NA
  #  for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))
    setwd(outtabledir)

    ct <- c("CC1", "CC2", "CC3", "CC4", "CC5")
    Rn = res[res$region=="cfanorth", 3:7]
    #Rn = as.matrix( res[ which(res$region=="cfanorth") , as.character(c(1:5)) ] )
    rownames(Rn) = years
    colnames(Rn) = ct
    print.xtable(Rn, type="latex", file="table.CC.large.north.obs.tex")
    HTML(Rn, file="table.CC.large.north.obs.html")

    #Rs = as.matrix( res[ which(res$region=="cfasouth") , as.character(c(1:5)) ] )
    Rs = res[res$region=="cfasouth", 3:7]
    rownames(Rs) = years
    colnames(Rs) = ct
    print.xtable(Rs, type="latex", file="table.CC.large.south.obs.tex")
    HTML(Rs, file="table.CC.large.south.obs.html")

    #Rx = as.matrix( res[ which(res$region=="cfa4x") , as.character(c(1:5)) ] )
    Rx = res[res$region=="cfa4x", 3:7]
    rownames(Rx) = years
    colnames(Rx) = ct
    print.xtable(Rx, type="latex", file="table.CC.large.4x.obs.tex")
    HTML(Rx, file="table.CC.large.4x.obs.html")

# ----------------------------------------
#  Percent soft from observed data  
  
    odb = odb0
    odb = odb[ which( odb$cw > 95 & odb$cw < 170 & odb$prodcd_id=="0" ) ,]  # commercial crab 
    years = sort( unique( odb$fishyr ) )


    res = NULL
    for (r in p$regions) {
    for (y in years) {
      out = proportion.soft (odb, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

    cnames = c("region", "fishyr", "pr.soft", "nsoft", "ntot")
    colnames(res) = cnames
    print(res)
    res = as.data.frame(res)

    for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))

    Rn = as.matrix( res[ which(res$region=="cfanorth") , as.character(c(1:3)) ] )
    rownames(Rn) = years
    latex(Rn, file="", title="", label="table.proportion.soft.north.obs", rowlabel="Year", cgroup="", na.blank=T, caption="Percent soft (by number) over time for N-ENS from at-sea-observed data.")

    Rs = as.matrix( res[ which(res$region=="cfasouth" ), as.character(c(1:5)) ] )
    rownames(Rs) = years
    latex(Rs, file="", title="", label="table.proportion.soft.south.obs", rowlabel="Year", cgroup="", na.blank=T, caption="Percent soft (by number) over time for S-ENS from at-sea-observed data.")

    Rx = as.matrix( res[ which(res$region=="cfa4x") , as.character(c(1:5)) ] )
    rownames(Rx) = years
    latex(Rx, file="", title="", label="table.proportion.soft.4x.obs", rowlabel="Year", cgroup="", na.blank=T, caption="Percent soft (by number) over time for CFA 4X from at-sea-observed data.")



# instars of interest: 11 and 12

# growth increment (assumming average weight in the midpoint of each increment)
  growth.11.to.12 =  predict.mass.g.from.CW.mm( mean(CW.interval.male(12)) ) - predict.mass.g.from.CW.mm (mean(CW.interval.male(11)) )
 
  # = 419 g
#  12to13 = ~450



  # Table of proportion discarded
    odb = observer.db("odb")
    regions = c("cfanorth", "cfasouth", "cfa4x")
    years = sort( unique( odb$fishyr ) )
    out = NULL
    for (r in regions) {
    for (y in years) {
      res = proportion.legal (odb, region=r, year=y)
      out = rbind(out, cbind( r, y, res[1], res[2], res[3] ) ) 
    } }
    out



# ---------------------------------------- USED 
#  Carapace condition from trawl data  >= 95mm CW  ... not kriged .. simple proportions
  
    det0 = snowcrab.db( DS="det.georeferenced" )
    det0$fishyr = det0$yr  ## the counting routine expectes this variable

    det = det0[ which( det0$cw >= 95 ) ,]  # commerical sized crab only
    years = sort( unique( det$yr ) )

    res = NULL
    for (r in p$regions) {
    for (y in years) {
      out = proportion.cc (det, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

    cnames = c("region", "fishyr", c(1:5), "ntot")
    colnames(res) = cnames
    print(res)
    res = as.data.frame(res)

    for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))
    (res)
















  # ------------------
  # counts of stations in each area

    # check towquality .. this should always == 1
    set = snowcrab.db("set.complete")
    if (length( unique( set$towquality) ) != 1 ) print("error -- not good tows")

    out = data.frame(yr=sort( unique(set$yr )) )
    for (reg in c("cfaall", "cfanorth", "cfasouth","cfa4x"  ) ) {
      d = filter.region.polygon(set[,c("lon","lat")], reg)
      e = as.data.frame( xtabs(~yr, data=set[d,])  )
      names(e) = c("yr", reg)
      e$yr = as.numeric(as.character(e$yr) )
      out = merge(out, e, by="yr", all=T)
    }
    print(out)

    x11()
    year = p$current.assessment.year
    setdata = set[ which(set$yr==year),]
    N = filter.region.polygon(setdata[,c("lon","lat")], "cfanorth")
    S = filter.region.polygon(setdata[,c("lon","lat")], "cfasouth")
    X = filter.region.polygon(setdata[,c("lon","lat")], "cfa4x")
    plot(setdata$lon, setdata$lat)
    points(setdata$lon[N], setdata$lat[N],col="red",pch=20)
    points(setdata$lon[S], setdata$lat[S],col="blue",pch=20)
    points(setdata$lon[X], setdata$lat[X],col="black",pch=20)




# % mat calculations:
#as above but

      loc = file.path(sc.R, "size.data")
      dir.create(path=loc, recursive=T, showWarnings=F)
      outfilename = paste( c("mi", "mm", "fi", "fm"), "rdata", sep=".")
      outfile = file.path(loc, paste(outfilename))
      for (f in  outfile) load(f)


           f.i = f.imm[which( rownames(f.imm)%in% sids ) ,]
           f.i.means = apply(X=f.i, MARGIN=2, FUN=mean)
           f.m = f.mat[which( rownames(f.mat)%in% sids ) ,]
           f.m.means = apply(X=f.m, MARGIN=2, FUN=mean)

           toplot = rbind(f.m.means, f.i.means)

 ii = as.data.frame(t(toplot))
 ii$cw = as.numeric(rownames(ii))
 ii$pmat = ii[,1]/ (ii[,1]+ii[,2]) * 100

plot(ii$cw, ii$pmat)
abline(h=50)

 str(ii)
#`data.frame':   70 obs. of  4 variables:
# $ f.m.means: num  0 0 0 0 0 ...
# $ f.i.means: num   2.80  6.19 20.05 24.29 74.11 ...
# $ cw       : num  12 14 16 18 20 22 24 26 28 30 ...
# $ pmat     : num  0 0 0 0 0 ...




# ----------------------------------------   NOT USED ____________
#  Carapace condition from trawl data  < 95mm CW  ... not kriged .. simple proportions
  
    det0 = snowcrab.db( DS="det.georeferenced" )
    det0$fishyr = det0$yr  ## the counting routine expectes this variable

    det = det0[ which( det0$cw < 95 ) ,]  # commerical sized crab only
    years = sort( unique( det$yr ) )

    res = NULL
    for (r in p$regions) {
    for (y in years) {
      out = proportion.cc (det, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

    cnames = c("region", "fishyr", c(1:5), "ntot")
    colnames(res) = cnames
    print(res)
    res = as.data.frame(res)

    for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))
    (res)





