
observer.load.rawdata = function(file) {
  out = NULL
  rawdata=NULL
  header = scan( file, what="", sep="\n", strip.white=F )  # there are control characters at the end of a file on occasion
  header = header[which(nchar(header)>140)]
  tmpfile =  file.path( tempdir(), make.random.string(".tmp.observer.csv"))

  write.table(header, file=tmpfile, append=F, row.names=F, col.names=F, quote=F )

  if (length(header) > 1) {
  }

  print("debug", 3)
  file.remove(tmpfile)

  return( rawdata )
}


