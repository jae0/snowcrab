
  tagging.example.usage = function() {
      
    loadfunctions( "snowcrab" )

    marked.file = "tags.1996_2001.csv"
    recaps.file = "recaptures.csv"

    recaps = get.recaps () #  alternate:(DS="raw")
    marked = get.marked () #  alternate:(DS="raw")
    marked2 =  read.table( file.path(project.datadirectory("snowcrab"), "data", "tagging", "tags_summary1993_2005.csv"), sep=";", header=T, as.is=T)

    move = get.move () #  alternate:(DS="redo")
      

    tmp0 = move[, c("lon0", "lat0")]
    names(tmp0) = c("lon", "lat")
    ss0 = filter.region.polygon(tmp0, region=region)

    tmp1 = move[, c("lon1", "lat1")]
    names(tmp1) = c("lon", "lat")
    ss1 = filter.region.polygon(tmp1, region=region)

    ss = unique(c(ss0, ss1)) 
    move = move[ ss, ]

  }


