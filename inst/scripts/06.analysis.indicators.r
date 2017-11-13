
  
  # analysis of "ecosystem indicators"
  
  if (!exists("year.assessment")) year.assessment=lubridate::year(Sys.Date())
  p = bio.snowcrab::load.environment( year.assessment=year.assessment )



	p$libs = unique( c( p$libs, project.library( 
    "stm", "parallel", "sorted.ordination", "emaf") ) )

  # setwd( project.datadirectory("emaf") )

  # not all are fully refreshed automatically .. they are just place holders for now

      groundfish = indicators.ts.db( db="groundfish.timeseries.redo" )
      snowcrab = indicators.ts.db( db="snowcrab.timeseries.redo")
      climate = indicators.ts.db (db="climate.redo" )
      shrimp = indicators.ts.db( db="shrimp.timeseries.redo")

      sar = indicators.ts.db( db="species.area.redo" )
      nss = indicators.ts.db( db="sizespectrum.redo" )
      metab = indicators.ts.db( db="metabolism.redo" )
      sc = indicators.ts.db( db="species.composition.redo" )

      economics = indicators.ts.db( db="economic.data.redo" )

      # hand constructed and updated .. TODO :: find solutions!
      #plankton = indicators.ts.db( db="plankton.timeseries.redo" )
      human = indicators.ts.db( db="demographics.redo" )


      #seals = indicators.ts.db( db="seal.timeseries.redo" )
      landedvalue = indicators.ts.db( db="landedvalue.annual.redo", ref.year=2008 )
      landings = indicators.ts.db( db="landings.annual.redo" )





# require( xlsReadWrite )
# data = read.xls( "mydata.xls", sheet="Sheet1" )
#
# for ( y in
# http://www.gov.ns.ca/finance/communitycounts/export.asp?bExcel=1&page=table_d17&dgroup=&dgroup1=&dgroup2=&dgroup3=&dgroup4=&yearid=2011&gnum=pro9012&gname=Nova%20Scotia&range=
#
# require( XLConnect )
# fn = "~/Downloads/estimates.xls"
# wb <- loadWorkbook( fn)
# data <- readWorksheet(wb)




  indic = indicators.ts.db( db="emaf.all.glue" )  # glue all time-series together
  # indic = indicators.ts.db( db="emaf.all" ) # load the glued version


 # indic$data$Nwells.drilled = cumsum.jae(indic$data$Nwells.drilled)
 # indic$data$seismic.2D = cumsum.jae(indic$data$seismic.2D)
 # indic$data$seismic.3D = cumsum.jae(indic$data$seismic.3D)


  t0 = 1960
  t1 = 2015

  # ordination of selected key factors
  indic = indicators.ts.db( db="emaf.all" )

  d = emaf.subset ( indic, type="keyfactors" )
#  save( d, file="/home/adam/tmp/ordin.rdata", compress=TRUE )

  Y = pca.analyse.data(d, t0, t1, fname=file.path(project.datadirectory("emaf"), "keyfactors" ) )


  sub = indic$data[, c("T_bottom_misaine", "SST_halifax", "ice_coverage.km.2", "Gulf.Stream.front.Ref.62lon", "T_sable_annual", "snowcrab.bottom.habitat.area", "snowcrab.kriged.R0.mass", "snowcrab.fishery.landings", "snowcrab.fishery.cpue", "groundfish.stratified.mean.temp" )]

  write.table(sub, file=file.path( project.datadirectory( "bio.snowcrab"), "research", "environ.management", "data.csv"), sep=";")


inn = names (indic$data)
for (i in .keyfactors) {
  if ( i %in% inn ) next()
  print (i)
}


  ## smaller subsets


  # human
  .human = c(indic$landings.totals.NS, indic$landedvalue.totals.NS, indic$human )
  .human = setdiff( .human, "No.Fish.processors" ) # this has no data yet
  Y = pca.analyse.data( indic$data, .human, t0, t1, fname=file.path(project.datadirectory("emaf"), "human") )

  # fishery landings
  .fishery = c(indic$landings.NS, indic$landings.totals.NS )
  Y = pca.analyse.data( indic$data, .fishery, t0, t1, fname=file.path(project.datadirectory("emaf"), "fishery" ))

  # fishery landed value
  .fishery.value = c(indic$landedvalue.NS, indic$landedvalue.totals.NS )
  Y = pca.analyse.data( indic$data, .fishery.value, t0, t1, fname=file.path(project.datadirectory("emaf"), "landedvalue" ))

  # fishery -- overall
  .fishery = c(indic$landedvalue.NS, indic$landedvalue.totals.NS, indic$landings.NS, indic$landings.totals.NS )
  Y = pca.analyse.data( indic$data, .fishery, t0=1970, t1, fname=file.path(project.datadirectory("emaf"), "fishery.overall" ))

  # snowcrab
  .snowcrab = c(indic$snowcrab, "groundfish.stratified.mean.totwgt.snowcrab", "groundfish.stratified.mean.totno.snowcrab" )
  Y = pca.analyse.data(indic$data, .snowcrab, t0, t1, fname=file.path(project.datadirectory("emaf"), "snowcrab" ))

  # climate
  .climate = c( indic$climate )
  Y = pca.analyse.data(indic$data, .climate, t0, t1, fname=file.path(project.datadirectory("emaf"), "climate" ))

  # ecosystem
  .ecosystem = c( indic$ecosystem )
  Y = pca.analyse.data(indic$data, .ecosystem, t0, t1, fname=file.path(project.datadirectory("emaf"), "ecosystem" ))








## ---------------
##  Todo : bayesian PCA

require(bmvm)
loadfunctions("bmvm")  ## mostly complete


 ...  incomplete

  # figures and tables related to fishery indices
  
  if (!exists("year.assessment")) year.assessment=lubridate::year(Sys.Date())
  p = bio.snowcrab::load.environment( year.assessment=year.assessment )



  project.library( "sorted.ordination")

  set = snowcrab.db("set.with.cat")
  allvars = c(
      "totno.all", "totno.male", "totno.male.com", "totno.male.ncom", "totmass.male.imm", "totmass.male.mat",
      "totmass.all", "totmass.male", "totmass.male.com", "totmass.male.ncom", "totno.male.imm", "totno.male.mat",

      "totno.male.soft", "totno.male.hard",
      "totmass.male.soft", "totmass.male.hard",

      "totno.female.soft", "totno.female.hard",
      "totmass.female.soft", "totmass.female.hard",

      "R0.no", "R1.no", "R2.no", "R3.no", "R4.no", "R5p.no", "dwarf.no",
      "R0.mass", "R1.mass", "R2.mass", "R3.mass", "R4.mass", "R5p.mass", "dwarf.mass",

      "totmass.female", "totmass.female.mat", "totmass.female.imm",
      "totmass.female.berried", "totmass.female.primiparous", "totmass.female.multiparous",

      "totno.female", "totno.female.mat", "totno.female.imm",
      "totno.female.berried", "totno.female.primiparous", "totno.female.multiparous",

      "totmass.male.com.CC1to2", "totmass.male.com.CC3to4", "totmass.male.com.CC5",
      "totno.male.com.CC1to2", "totno.male.com.CC3to4", "totno.male.com.CC5",

      "cw.fem.mat.mean", "cw.fem.imm.mean",
      "cw.male.mat.mean", "cw.male.imm.mean",

      "cw.mean", "cw.comm.mean", "cw.notcomm.mean",

      "sexratio.imm", "sexratio.mat", "sexratio.all",

      "z", "t",

      "uniqueid"
    )

    x = set[, allvars]
    for (v in allvars)  x[,v] = recode.variable(x[,v], v)$x  # normalise the data where required

    rownames(x) = x$uniqueid
    x$uniqueid = NULL
    p = pca.modified(x)




