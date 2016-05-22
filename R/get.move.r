 
  get.move = function( redo=F ) {
    move = NULL
    require(chron) 
    tags.datadir= file.path( project.datadirectory("snowcrab"), "data", "tagging" )
    
    fn =file.path( tags.datadir, paste("move", "rdata", sep="." ) )

    if ( !redo ) {
      load( fn )
      return (move) 
    }
    
    recaps = get.recaps (DS="file") 
    marked = get.marked (DS="file")
    marked2 =  read.table( file.path(tags.datadir, "tags_summary1993_2005.csv" ), sep=";", header=T, as.is=T)

    v0 = c("tagID", "jul", "lon", "lat", "cw", "ch", "cc", "z.fm", "area" )
    v1 = c("tagID", "jul", "lon", "lat", "cw", "ch", "cc", "fisherman", "z.fm", "duro", "Comments")
    move = merge(x=recaps[,v1], y=marked[,v0], by="tagID" , sort=F, all.x=T, all.y=F, suffixes=c("1","0") )
    
    strange = NULL
    strange = which(move$lon1 < -66)
    strange = c(strange, which(move$lon1 < -62 & move$lat1 > 45 & move$lat0 < 45) )
      
    move = move[ -strange, ]
    dim(move)

    # fix data with no starting positions   
    no.data = which( is.na(move$lon0 + move$lat0) )
    for (i in no.data) {
      studyid = determine.origin(move[i,"tagID"])
      if (!is.null(studyid)){
          s = marked2[studyid, c("lon", "lat", "yr")]
          move$lon0[i] = s$lon
          move$lat0[i] = s$lat
          if ( is.na(move$jul0[i])) move$jul0[i] = chron( dates.=paste(s$yr, 12, 1, sep="-"),
                            format=c(dates="y-m-d") )
        }
    }
    move$dt = move$jul1 - move$jul0
    r = which(move$dt > -365 & move$dt < 0)
    move[r, "dt"] = 4*31  # caught in same season .. assume length of fishing season ~ 4 months

    move$dx = great.circle.distance.vector.vector(
      move[, c("lon0", "lat0")], move[, c("lon1", "lat1")], R=6367.436 
    )
   
    save(move,  file=fn, compress=T)
    return (move)
  }


