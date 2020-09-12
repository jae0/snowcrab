
  # Figures and tables obtained after completion of data assimilation and processing up to the end of "01.snowcrab.r"
require(aegis)


year.assessment = 2019

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
# loadfunctions('bio.snowcrab')


  # --- FIGURES ----

  # ------------------------------------------
   # Time-series: Fisheries landings
  #BZ TODO get these to x the anuary meetings script exactly as that lines up with Improptu reports
  #for now, the plots produced are accurate.
   figure.landings.timeseries( yearmax=p$year.assessment, outdir=file.path( p$annual.results,  "timeseries","fishery"), outfile="landings.ts", outfile2="landings.ts.sm" )

  # ------------------------------------------
  # Time-series: Fisheries effort
   figure.effort.timeseries( yearmax=p$year.assessment, outdir=file.path( p$annual.results,"timeseries", "fishery"), outfile="effort.ts", outfile2="effort.ts.sm" )

  # ------------------------------------------
  # Time-series: Fisheries CPUE
   figure.cpue.timeseries( yearmax=p$year.assessment, outdir=file.path( p$annual.results,"timeseries", "fishery"), outfile="cpue.ts", outfile2="cpue.sm.ts" )

  # ------------------------------------------
  # Size frequency distributions, broken down by moult category from at-sea observed data

    figure.observed.size.freq( regions = c("cfanorth", "cfasouth", "cfa4x"), years="all", outdir=file.path( p$annual.results, "figures", "size.freq", "observer")  )

  # ------------------------------------------
  # Size-frequency distributions of snow crab cw from trawl data, broken down by maturity classes
    histograms.size.maturity.update( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  redo.data=T )
    histograms.size.maturity.single.area( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  area='cfa4x',redo.data=T ) #area = cfanorth, cfasouth of cfa4x


  # ------------------------------------------
  # Timeseries of all survey variables
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), variables="R0.mass", plotyears=2004:p$year.assessment) # just R0 to see
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"),variables=c("sexratio.all","sexratio.mat","sexratio.imm"))
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"),plotyears=2004:p$year.assessment) # all variables
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "observer"),plotyears=2004:p$year.assessment,type='observer')

   figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"),plotyears=2004:p$year.assessment,type='groundfish.t') # groundfish survey temperature
  #-----------------------------------------------

  #Timeseries: geometric mean biomass of by-catch from snow crab survey

  # predators and competitors
    #cod, haddock, halibut, plaice, wolfish, thornyskate, smoothskate, winterskate, northernshrimp, jonahcrab, lessertoadcrab
  species = c(10, 11, 30, 40, 201, 50, 2521, 2511, 202, 204, 2211)
  figure.timeseries.bycatch(p=p, species=species, plotyears=2004:p$year.assessment, outdir=file.path(p$annual.results,"timeseries", "survey"))


  # ------------------------------------------
  # Map:  Interpolated mean/geometric mean of various variables in the set data table
#  p$do.parallel=F
  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )

  outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )


  # just for the roadshow
    map.set.information( p=p, outdir=outdir, variables=c('totmass.male.com', 'totmass.female.mat'),mapyears=p$mapyears)
    map.set.information( p=p, variables='t',mapyears=p$mapyears,outdir=outdir,log.variable=F,add.zeros=F,theta=100)

    # bycatch (geometric means)
    bc.vars = c(paste("ms.mass",species,sep='.'),paste("ms.no",species,sep='.'))
    outdir.bc= file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )
    map.set.information( p, variables=bc.vars, mapyears=p$mapyears, outdir=outdir.bc,probs=c(0,0.975)) #


    # all variables (geometric means)
    #map.set.information( p, outdir=outdir) # takes a long time

    # Means
    # variables that shouldn't be logged
    set = snowcrab.db( p=p, DS="set.biologicals")
    variables = bio.snowcrab::snowcrab.variablelist("all.data")
    variables = intersect( variables, names(set) )

    nolog.variables = c("t","z","sexratio.all","sexratio.mat","sexratio.imm","julian","julian.compressed", variables[grep("cw",variables)])
    map.set.information( p=p, variables=nolog.variables,outdir=outdir,log.variable=F,add.zeros=F,theta=35)
    # logit transform for ratios
    map.set.information( p=p, variables=c("sexratio.all","sexratio.mat","sexratio.imm"),outdir=outdir,log.variable=F,add.zeros=F,theta=40)

    # Geometric Means
    # all except variables that shouldn't be logged
    mass.vars = variables[!variables%in%nolog.variables][grep('mass',variables[!variables%in%nolog.variables])]
    no.vars = variables[!variables%in%nolog.variables][grep('no',variables[!variables%in%nolog.variables])]
    map.set.information( p=p, variables= mass.vars,outdir=outdir)
    map.set.information( p=p, variables= no.vars,outdir=outdir,probs=c(0,0.975))



  # ------------------------------------------
  # Map: Survey locations

    map.survey.locations( p=p, basedir=file.path(p$project.outputdir, "maps", "survey.locations"),  newyear=F, map.method="lattice"  )
  #  map.survey.locations( p, basedir=file.path(p$project.outputdir, "maps", "survey.locations"),  newyear=F, map.method="googleearth"  )

  # ------------------------------------------
  # Map: Observer locations
    map.observer.locations( p=p, basedir=file.path(p$project.outputdir, "maps","observer.locations" ), newyear=F , map.method="lattice"  )

  # ------------------------------------------
  # Map: Logbook recorded locations
    map.logbook.locations( p=p, basedir=file.path(p$project.outputdir, "maps","logbook.locations" ), newyear=F , map.method="lattice"  )



  # ------------------------------------------
  # Map: Logbook data
  outdir = file.path( p$project.outputdir, "maps", "logbook","snowcrab","annual" )
  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )

  map.fisheries.data( p=p, variable= 'effort', outdir=outdir, FUN=sum, probs=c(0,0.975))
  map.fisheries.data( p=p, variable= 'cpue', outdir=outdir, FUN=mean, probs=c(0,0.975))
  map.fisheries.data( p=p, variable= 'landings', outdir=outdir, FUN=sum, probs=c(0,0.975))


  # --- TABLES ----
  # TODO-BZ add functionality for tables to be saved as pdf
  # add tab.4.tex.r function

  require(gridExtra)
  library("xtable")
  library("R2HTML")

  odb0 = observer.db("odb")
  regions = c("cfanorth", "cfasouth", "cfa4x")
  nregions = length(regions)

  #------------------------------------------------
  #Fisheries statistics per region
  tabledir = file.path(project.datadirectory("bio.snowcrab"), "data", "fisheries")
  outtabledir= file.path(project.datadirectory("bio.snowcrab"), "assessments", p$year.assessment, "tables", "logbook")
  if(!dir.exists(tabledir)) dir.create(tabledir, recursive =T)
  if(!dir.exists(outtabledir)) dir.create(outtabledir, recursive =T)

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

    outtabledir= file.path(project.datadirectory("bio.snowcrab"), "assessments", p$year.assessment, "tables", "observer")

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

    det0 = snowcrab.db( p=p, DS="det.georeferenced" )
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

    HTML(res, file="table.CC.large.survey.html")

  # ------------------
  # counts of stations in each area

    # check towquality .. this should always == 1
    set = snowcrab.db(p=p, DS="set.clean")
    if (length( unique( set$towquality) ) != 1 ) print("error -- not good tows")

    out = data.frame(yr=sort( unique(set$yr )) )
    for (reg in c("cfaall", "cfanorth", "cfasouth","cfa4x"  ) ) {
      d = polygon_inside(set[,c("lon","lat")], reg)
      e = as.data.frame( xtabs(~yr, data=set[d,])  )
      names(e) = c("yr", reg)
      e$yr = as.numeric(as.character(e$yr) )
      out = merge(out, e, by="yr", all=T)
    }
    print(out)

    plot.new()
    year = p$year.assessment
    setdata = set[ which(set$yr==year),]
    N = polygon_inside(setdata[,c("lon","lat")], "cfanorth")
    S = polygon_inside(setdata[,c("lon","lat")], "cfasouth")
    X = polygon_inside(setdata[,c("lon","lat")], "cfa4x")
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

    det0 = snowcrab.db( p=p, DS="det.georeferenced" )
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





  ##########################################################################
  ###############################  Retired figures #########################



    # ------------------------------------------
  # Map: Scotian Shelf with CFA lines and labels  .. using gmt
  # this is the basemap from map.r which is then post-labelled in sodipodi
  #  p$outdir = file.path(p$annual.results,"figures")
  #  p$outfile.basename = file.path(p$outdir, "map.CFAs")
  #  map.basemap.with.cfa.lines( p, conversions=c("ps2png")  )


  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
   # habitat.usage( usevar="totno.all", covariate="depth", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    #habitat.usage( usevar="totno.all", covariate="temperature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
   # habitat.usage( usevar="totno.all", covariate="bottom.slope", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    #habitat.usage( usevar="totno.all", covariate="bottom.curvature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    #habitat.usage( usevar="totno.all", covariate="substrate", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Timeseries: Larval brachyura from the SSIP data
    ##figure.timeseries.larvae( outdir=file.path(p$project.outputdir, "timeseries", "larvae") )

  # ------------------------------------------
  # Growth as a a function of instar for Scotian Shelf snow crab
    figure.growth.instar( outdir=file.path(p$project.outputdir, "growth") )


  # ------------------------------------------
  # Map: Larval distributions from the Scotian Shelf Ichtyoplankton Program data
    map.larvae( p=p, outdir=file.path(p$project.outputdir, "maps", "larvae"), conversions=conversions )


    # ------------------------------------------
    # Map: Spatial representation of maturity patterns of snow crab
    #MG Not sure we use these maps either, check with Adam and Jae
    # map.maturity( p, outdir=file.path(p$project.outputdir, "maps", "maturity"), newyear=T )
