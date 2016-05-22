 map.survey.locations = function(p, basedir, newyear=T, map.method="GMT" ) { 


    set = snowcrab.db( DS="set.clean")
    years = sort( unique( set$yr ) )
    if (newyear) years = p$current.assessment.year

    if ( map.method=="GMT" ) {
      
      set = set[, c("yr", "lon", "lat")]
      set = set[ is.finite( rowSums(set) ) ,]
 
      p$psxyoptions = "-Sc0.1c -G20"  # Sc = circle with size 0.1cm, G is color/grayscale
      p$basedir = basedir
     
      for (y in years) {
        toplot = set[ which(set$yr==y), c("lon", "lat")]
        p$outfile.basename = file.path(p$basedir, paste("survey.locations", y, sep=".") )
        print(  p$outfile.basename )
        gmt.xyplot ( p, toplot, y, conversions=p$conversions  )
      }

      pause(30)
      files.to.delete = list.files( p$basedir, "^survey.locations.*.ps$", all.files=T, full.names=T, recursive=T) 
      remove.files ( files.to.delete ) 
    
    }
    
    if (map.method=="lattice" ) {
      
      set = set[, c("yr", "plon", "plat")]
      set = set[ is.finite( rowSums(set) ) ,]

      for (y in years) {
        toplot = set[ which(set$yr==y), c("plon", "plat")]
        annot = paste ("Survey locations", y)
        fn = paste("survey.locations", y, sep=".") 
        print(fn)
        map( toplot, cfa.regions=T, depthcontours=T, annot=annot, fn=fn, loc=basedir, corners=p$planar.corners )
       
      }
    }

    if (map.method=="googleearth" ) {
      
      loadfunctions( "plottingmethods" )
  
      outfile =  file.path( basedir, "survey.locations.googleearth.kml" ) 
      set$desc = paste( 
              "Trip:", set$trip, "\n", 
              "Set:", set$set, "\n", 
              "Station:", set$station, "\n",
              "Depth:", set$z, "\n",
              "Temperature:", set$t
              ) 

      coords = c("lon", "lat", "elevation")
      set$elevation = 0
#      set = set[ is.finite( rowSums(set[,coords]) ) ,]

      # start kml document
      con = kml.start( outfile, "Survey stations" )
      
      # define point styles/colours, etc
      kml.placemark.make( con, item="style", style.id="pin.red", colour="c0ffffff", scale=0.30, 
        href='http://nssnowcrab.googlepages.com/reddot.png' )  # red dot

      # main folder start
      kml.folder.start( con, folder.name="Scotian Shelf snow crab survey stations", 
        desc="Asssessment tools for the Scotian Shelf snow crab (Bedford Institute of Oceanography)" 
      )
        for (y in years) {
          toplot = set[ which(set$yr==y), ]
          kml.folder.start( con, y )
            for ( i in 1:nrow(toplot) ) { # low defined at top
              kml.placemark.make( con, desc=toplot[i,"desc"], style.id="pin.red", x=toplot[i, coords] ) 
            }
          kml.folder.end( con )
        
        }

      kml.folder.end( con )
      kml.end( con )


      print( "Google earth file has been completed and saved as the following: ")
      print( outfile )
  
    }

    return ("Done" )
  }


