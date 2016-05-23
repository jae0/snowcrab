 

  require(chron)

	loadfunctions( c( "spacetime", "utility", "parallel", "sorted.ordination", "indicators") ) 
 
  setwd( project.datadirectory("indicators") )
        
  # not all are fully refreshed automatically .. they are just place holders for now
      
      groundfish = indicators.db( db="groundfish.timeseries.redo" )
      snowcrab = indicators.db( db="snowcrab.timeseries.redo") 
      climate = indicators.db (db="climate.redo" )
      shrimp = indicators.db( db="shrimp.timeseries.redo")

      sar = indicators.db( db="species.area.redo" )
      nss = indicators.db( db="sizespectrum.redo" )
      metab = indicators.db( db="metabolism.redo" )
      sc = indicators.db( db="species.composition.redo" )

      economics = indicators.db( db="economic.data.redo" )

      # hand constructed and updated .. TODO :: find solutions!
      #plankton = indicators.db( db="plankton.timeseries.redo" )
      human = indicators.db( db="demographics.redo" )
     
      
      #seals = indicators.db( db="seal.timeseries.redo" ) 
      landedvalue = indicators.db( db="landedvalue.annual.redo", ref.year=2008 )
      landings = indicators.db( db="landings.annual.redo" )
      

  
  
  
  
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




  indic = indicators.db( db="indicators.all.glue" )  # glue all time-series together
  # indic = indicators.db( db="indicators.all" ) # load the glued version
  

 # indic$data$Nwells.drilled = cumsum.jae(indic$data$Nwells.drilled)
 # indic$data$seismic.2D = cumsum.jae(indic$data$seismic.2D)
 # indic$data$seismic.3D = cumsum.jae(indic$data$seismic.3D)


  t0 = 1960
  t1 = 2015

  # ordination of selected key factors
  indic = indicators.db( db="indicators.all" )  
  
  d = indicators.subset ( indic, type="keyfactors" )
#  save( d, file="/home/adam/tmp/ordin.rdata", compress=TRUE )
  
  Y = pca.analyse.data(d, t0, t1, fname=file.path(project.datadirectory("indicators"), "keyfactors" ) )
  

  sub = indic$data[, c("T_bottom_misaine", "SST_halifax", "ice_coverage.km.2", "Gulf.Stream.front.Ref.62lon", "T_sable_annual", "snowcrab.bottom.habitat.area", "snowcrab.kriged.R0.mass", "snowcrab.fishery.landings", "snowcrab.fishery.cpue", "groundfish.stratified.mean.temp" )]

  write.table(sub, file=file.path( project.datadirectory( "snowcrab"), "research", "environ.management", "data.csv"), sep=";")


inn = names (indic$data) 
for (i in .keyfactors) {
  if ( i %in% inn ) next()
  print (i)
}


  ## smaller subsets


  # human 
  .human = c(indic$landings.totals.NS, indic$landedvalue.totals.NS, indic$human ) 
  .human = setdiff( .human, "No.Fish.processors" ) # this has no data yet
  Y = pca.analyse.data( indic$data, .human, t0, t1, fname=file.path(project.datadirectory("indicators"), "human") )

  # fishery landings
  .fishery = c(indic$landings.NS, indic$landings.totals.NS ) 
  Y = pca.analyse.data( indic$data, .fishery, t0, t1, fname=file.path(project.datadirectory("indicators"), "fishery" ))
  
  # fishery landed value
  .fishery.value = c(indic$landedvalue.NS, indic$landedvalue.totals.NS ) 
  Y = pca.analyse.data( indic$data, .fishery.value, t0, t1, fname=file.path(project.datadirectory("indicators"), "landedvalue" ))

  # fishery -- overall
  .fishery = c(indic$landedvalue.NS, indic$landedvalue.totals.NS, indic$landings.NS, indic$landings.totals.NS ) 
  Y = pca.analyse.data( indic$data, .fishery, t0=1970, t1, fname=file.path(project.datadirectory("indicators"), "fishery.overall" ))

  # snowcrab 
  .snowcrab = c(indic$snowcrab, "groundfish.stratified.mean.totwgt.snowcrab", "groundfish.stratified.mean.totno.snowcrab" ) 
  Y = pca.analyse.data(indic$data, .snowcrab, t0, t1, fname=file.path(project.datadirectory("indicators"), "snowcrab" ))

  # climate
  .climate = c( indic$climate )
  Y = pca.analyse.data(indic$data, .climate, t0, t1, fname=file.path(project.datadirectory("indicators"), "climate" ))

  # ecosystem
  .ecosystem = c( indic$ecosystem )
  Y = pca.analyse.data(indic$data, .ecosystem, t0, t1, fname=file.path(project.datadirectory("indicators"), "ecosystem" ))




