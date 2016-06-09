
  p = bio.snowcrab::initialise.local.environment( current.assessment.year=2016)

	p$libs = unique( c( p$libs, bioLibrary( "spacetime", "bio.utilities", "parallel", "sorted.ordination", "bio.indicators") ) )

  # setwd( project.datadirectory("bio.indicators") )

  # not all are fully refreshed automatically .. they are just place holders for now

      groundfish = bio.indicators.db( db="groundfish.timeseries.redo" )
      bio.snowcrab = bio.indicators.db( db="bio.snowcrab.timeseries.redo")
      climate = bio.indicators.db (db="climate.redo" )
      shrimp = bio.indicators.db( db="shrimp.timeseries.redo")

      sar = bio.indicators.db( db="species.area.redo" )
      nss = bio.indicators.db( db="sizespectrum.redo" )
      metab = bio.indicators.db( db="metabolism.redo" )
      sc = bio.indicators.db( db="species.composition.redo" )

      economics = bio.indicators.db( db="economic.data.redo" )

      # hand constructed and updated .. TODO :: find solutions!
      #plankton = bio.indicators.db( db="plankton.timeseries.redo" )
      human = bio.indicators.db( db="demographics.redo" )


      #seals = bio.indicators.db( db="seal.timeseries.redo" )
      landedvalue = bio.indicators.db( db="landedvalue.annual.redo", ref.year=2008 )
      landings = bio.indicators.db( db="landings.annual.redo" )






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




  indic = bio.indicators.db( db="bio.indicators.all.glue" )  # glue all time-series together
  # indic = bio.indicators.db( db="bio.indicators.all" ) # load the glued version


 # indic$data$Nwells.drilled = cumsum.jae(indic$data$Nwells.drilled)
 # indic$data$seismic.2D = cumsum.jae(indic$data$seismic.2D)
 # indic$data$seismic.3D = cumsum.jae(indic$data$seismic.3D)


  t0 = 1960
  t1 = 2015

  # ordination of selected key factors
  indic = bio.indicators.db( db="bio.indicators.all" )

  d = bio.indicators.subset ( indic, type="keyfactors" )
#  save( d, file="/home/adam/tmp/ordin.rdata", compress=TRUE )

  Y = pca.analyse.data(d, t0, t1, fname=file.path(project.datadirectory("bio.indicators"), "keyfactors" ) )


  sub = indic$data[, c("T_bottom_misaine", "SST_halifax", "ice_coverage.km.2", "Gulf.Stream.front.Ref.62lon", "T_sable_annual", "bio.snowcrab.bottom.habitat.area", "bio.snowcrab.kriged.R0.mass", "bio.snowcrab.fishery.landings", "bio.snowcrab.fishery.cpue", "groundfish.stratified.mean.temp" )]

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
  Y = pca.analyse.data( indic$data, .human, t0, t1, fname=file.path(project.datadirectory("bio.indicators"), "human") )

  # fishery landings
  .fishery = c(indic$landings.NS, indic$landings.totals.NS )
  Y = pca.analyse.data( indic$data, .fishery, t0, t1, fname=file.path(project.datadirectory("bio.indicators"), "fishery" ))

  # fishery landed value
  .fishery.value = c(indic$landedvalue.NS, indic$landedvalue.totals.NS )
  Y = pca.analyse.data( indic$data, .fishery.value, t0, t1, fname=file.path(project.datadirectory("bio.indicators"), "landedvalue" ))

  # fishery -- overall
  .fishery = c(indic$landedvalue.NS, indic$landedvalue.totals.NS, indic$landings.NS, indic$landings.totals.NS )
  Y = pca.analyse.data( indic$data, .fishery, t0=1970, t1, fname=file.path(project.datadirectory("bio.indicators"), "fishery.overall" ))

  # bio.snowcrab
  .bio.snowcrab = c(indic$bio.snowcrab, "groundfish.stratified.mean.totwgt.bio.snowcrab", "groundfish.stratified.mean.totno.bio.snowcrab" )
  Y = pca.analyse.data(indic$data, .bio.snowcrab, t0, t1, fname=file.path(project.datadirectory("bio.indicators"), "bio.snowcrab" ))

  # climate
  .climate = c( indic$climate )
  Y = pca.analyse.data(indic$data, .climate, t0, t1, fname=file.path(project.datadirectory("bio.indicators"), "climate" ))

  # ecosystem
  .ecosystem = c( indic$ecosystem )
  Y = pca.analyse.data(indic$data, .ecosystem, t0, t1, fname=file.path(project.datadirectory("bio.indicators"), "ecosystem" ))




