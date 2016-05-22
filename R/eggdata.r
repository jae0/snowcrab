
eggdata = function( data.loc=".", redo.data=F, data.file ) {

  if (redo.data) {
    
    require( chron)

    # Data from Jae Jan 14
    load(file.path("det.rdata"))
    load(file.path("set_biologicals.rdata"))

    # Need to add time back into the SNCRABDETAILS data set
    set$sdate = set$t0     #'sdate' is 't0' in new data set
    set.with.time = set[ , c("trip", "set", "sdate", "t") ]  # trip and set act as 'markers' to align data
    det = merge( det, set.with.time, by=c("trip", "set"), all.x=T, all.y=F)

    sc=set
    sdt=det

    # There are NA's in the date, this fixes the problem
    xyr = (substring( sdt$trip, 6, 9))
    xmon =  (substring( sdt$trip, 4, 5))
    xday =  (substring( sdt$trip, 2, 3))
    date.string = paste( xyr, xmon, xday, sep="/")
     
    date.from.trip = chron( dates.=date.string, format=c( dates = "y/m/d" ), out.format=c(dates="year-m-d")  )
     
    i = which(!is.finite(sdt$sdate))
    sdt$chrono = as.chron(sdt$sdate)
    sdt$chrono[i] = date.from.trip[i]
     
    sdt$yr = years(sdt$chrono)
    sdt$month = months(sdt$chrono)
    sdt$julianday =  convert.datecodes(sdt$chrono, "julian")
    sdt$dummy = 1

    sc$chronosc = as.chron(sc$sdate)
    sc$yr = years(sc$chronosc)
    sc$month = months(sc$chronosc)
    sc$day = days(sc$chronosc)
    sc$dummy = 1     # This is to use with 'xtabs'
    sc$julianday = convert.datecodes(sc$chronosc, "julian")
    sdt=sdt[!is.na(sdt$sex),]   # Takes all good data, i.e. that which is not NA
    sdt=sdt[!is.na(sdt$cw),]

    save( sc, file= file.path(data.loc, "sc.rdata"), compress=T )
    save( sdt, file =file.path(data.loc, "sdt.rdata"), compress=T )
    rm(set,det)
    return ("Files saved to work directory")
  }
  
  if (data.file =="sc") {
    load( file.path(data.loc, "sc.rdata") )
    return( sc )
  }

  if( data.file== "sdt" ) {
    load( file.path(data.loc, "sdt.rdata") )
    return(sdt)
  }

}


