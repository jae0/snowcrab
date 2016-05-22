
  get.marked = function(DS="file" ) { 
    
    tags.datadir= file.path( project.datadirectory("snowcrab"), "data", "tagging" )
    marked.file="tags.1996_2001.csv"
    marked = NULL
    
    if (DS=="redo") {
      marked = read.table( file.path( tags.datadir, marked.file), sep=";", header=T, as.is=T)
      marked$Ncrabs = NULL
      f = which(marked$lon>-55)
      marked[f,] = NA
      marked$jul = chron::chron(dates.=marked$date, format=c(dates="m/d/y"), out.format=c(dates="year-m-d") )
      save(marked, file=file.path(tags.datadir, "marked.Rdata"), compress=T)
    }
    if (DS=="file") load( file=file.path(tags.datadir, "marked.Rdata"))
    
    return (marked)

  }
  

