 
 
  get.recaps = function(DS="file") {
    
    tags.datadir= file.path( project.datadirectory("snowcrab"), "data", "tagging" )

    recaps.file="recaptures.csv"
    recaps = NULL
    
    if (DS=="redo") {
      recaps = read.table( file.path(tags.datadir, recaps.file), sep=";", header=T, as.is=T)
      f = which(recaps$lon>-50 & recaps$lon < 100 ) 
      recaps[f,"lon"] = -recaps[f,"lon"] 
      recaps = recaps[!is.na(recaps$yr) ,]
      
      f = which(recaps$month>12)
      months = recaps[f, "month"]
      recaps[f, "month"] = recaps[f, "date"]
      recaps[f, "date"] = months
      
      f = which(is.na(recaps$month))
      recaps[f, "month"] = 12 # assume a December capture
     
      f = which(recaps$date>31)
      recaps[f, "date"] = 1  # assume first day of the month

      f = which(is.na(recaps$date))
      recaps[f, "date"] = 1  # assume first day of the month
      
      recaps$jul = lubridate::ymd( dates.=paste(recaps$yr, recaps$month, recaps$date, sep="-"))
      recaps$jul = chron::chron( dates.=paste(recaps$yr, recaps$month, recaps$date, sep="-"),
        format=c(dates="y-m-d"), out.format=c(dates="year-m-d") )

      cc3 = which(recaps$cc %in% c(
       "Intermediate", "Legal/hard(released)", "clean",
       "intermediate", "propre", " Legal (released)         ", 
       "Clean / propre", "Intermediate/interm\351diaire", 
       "Interm\351diate", "good condition", 
       "very good condition", "Hard shell", "Legal (released)",
       "inter"))
      ccM = which(recaps$cc %in% c( "dirty", "Dirty / sale", "3M", "3m", "Mossy/Dirty", "mousseux",
                                    "intermediate (3M)", "sale", "mossy"))
      ccD = which(recaps$cc %in% c( "Undersize(released)      ", "Undersize(released)       ",
                                    " Undersize(released)      " ))
      cc5 = which(recaps$cc %in% c( "4 et 5" ))
          
      recaps[cc3,"cc"] = "3"
      recaps[ccM,"cc"] = "M"
      recaps[ccD,"cc"] = "Dw"
      recaps[cc5,"cc"] = "5"
      recaps$tagID = gsub("[ -]", "", tolower(recaps$tagID))
      save(recaps, file=file.path(tags.datadir, "recaps.Rdata"), compress=T)
    }
    if (DS=="file") load(file=file.path(tags.datadir, "recaps.Rdata"))
   
    return(recaps)
  }


