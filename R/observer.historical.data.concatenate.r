
  observer.historical.data.concatenate = function(type, outfile) {

    rootdir = file.path( project.datadirectory ("snowcrab"), "data", "observer" )

    directories =  file.path( rootdir, "archive", c(2001:2003) )
    directories =  file.path( rootdir, "archive", c(1996, 1999:2000) )
    directories =  file.path( rootdir, "archive", c(1997, 1998) )

    fl = list.files(path=directories, pattern="[*.txt]$", full.names=T, recursive=T)

    if (type=="seasamples") {
      fl = fl[ grep( "sea sample", fl, ignore.case=T ) ]
    } else if (type=="portsamples") {
      fl = fl[ grep( "port sample", fl, ignore.case=T ) ]
    }

    temp.output = "tmp.dump"
    if (file.exists( temp.output ) )  file.remove(temp.output)
    for (i in fl) {
      print ( paste("reading",i) )
      out = NULL
      rawdata=NULL
      read0 = scan( file=i, what="", sep="\n", strip.white=F )  # there are control characters at the end of a file on occasion
      read0 = read0[which(nchar(read0)>140)]
      write.table( read0, file=temp.output, append=T, row.names=F, col.names=F, quote=F, sep=";" )
    }

    # data format is column-width structured
    dataformat = observer.data.formats(2003)
    dataformat = observer.data.formats(1997)

    rawdata = read.fwf( file=temp.output, widths=dataformat$w, as.is=T, colClasses="character", header=F, strip.white=F)
    colnames(rawdata) = dataformat$v
    rawdata = rawdata[, dataformat$keep]

    save( rawdata, file="seasamples.1996.1999.2000.rdata", compress=T)
    save( rawdata, file="portsamples.1996.1999.2000.rdata", compress=T)

    save( rawdata, file="seasamples.1997.1998.rdata", compress=T)
    save( rawdata, file="portsamples.1997.1998.rdata", compress=T)

    save( rawdata, file=outfile, compress=T)

    save( rawdata, file=outfile, compress=T)
    if (file.exists( temp.output ) )  file.remove(temp.output)

    return( "Done" )

  }


