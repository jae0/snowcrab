	
	loadfunctions( "plottingmethods" ) 
  loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  

  workdir = file.path( project.datadirectory("snowcrab"), "R" )
   
  redo.background.kml = F
  if ( redo.background.kml ) {
    # make snow crab background layers to a separate file
    backgroundfile = file.path( workdir, "snowcrab_layers.kml" ) 
    kml.snowcrab.background( backgroundfile )
  }

  outfile =  file.path( workdir, "mark.recaptures.kml" ) 

  ss = get.move ( )
    region = "cfaall"
    tmp0 = ss[, c("lon0", "lat0")]
      names(tmp0) = c("lon", "lat")
      ss0 = filter.region.polygon(tmp0, region=region)
    tmp1 = ss[, c("lon1", "lat1")]
      names(tmp1) = c("lon", "lat")
      ss1 = filter.region.polygon(tmp0, region=region)
  ss = ss[ unique( c( ss0,ss1) ) , ]
  
  # googleearth does not like '&' .. remove them
  ss$fisherman = gsub( '&', 'and',  ss$fisherman, perl = TRUE, useBytes=T ) 
  ss$Comments = gsub( '&', 'and',  ss$Comments, perl = TRUE, useBytes=T ) 

  ss$elevation0 = - 1.8288 * as.numeric(ss$z.fm0)
  ss$elevation1 = - 1.8288 * as.numeric(ss$z.fm1)
  ss$elevation0[ !is.finite( ss$elevation0) ] = 0
  ss$elevation1[ !is.finite( ss$elevation1) ] = 0

  ss$desc0 = paste( 
    "Date:", ss$jul0 ) 
  
  ss$desc1 = paste( 
    "Date:", ss$jul0, "\n", 
    "Fisher:", ss$fisherman, "\n", 
    "Comments:", ss$Comments ) 

  coords0 = c("lon0", "lat0", "elevation0") # starting 
  coords1 = c("lon1", "lat1", "elevation1") # ending
  
   # start kml document
  con = kml.start( outfile, "snow crab mark-recaptures" )
    
    # define point styles/colours, etc
    kml.placemark.make( con, item="style", style.id="pin.red", colour="c0ffffff", scale=0.2, 
      href='http://nssnowcrab.googlepages.com/reddot.png' )  # red dot
    kml.placemark.make( con, item="style", style.id="pin.yellow", colour="a0ffffff", scale=0.15, 
      href='http://nssnowcrab.googlepages.com/yellowdot.png' )  # yellow dot
    kml.line.make( con, item="style", style.id="line.movement", line.colour="88ffffff", line.width=2 )
      
     
    kml.folder.start( con, folder.name="Tagging locations" )
      tags = ss [- which( duplicated( ss[, coords0] )) ,] 
      tags = tags [ which( is.finite( rowSums(tags[,coords0] ))) ,] 
       
      for ( i in 1:nrow(tags) ) { # low defined at top
        kml.placemark.make( con, desc=tags[i,"desc0"], style.id="pin.red", x=tags[i, coords0] ) 
      }
    kml.folder.end (con)
   
    kml.folder.start( con, folder.name="Recapture locations" )
      recap = ss[ which( is.finite( rowSums(ss[, coords1] ))) ,] 
      for ( i in 1:nrow(recap) ) { # low defined at top
        kml.placemark.make( con, desc=recap[i,"desc1"], style.id="pin.yellow", x=recap[i, coords1] ) 
      }
    kml.folder.end (con)
     
    kml.folder.start( con, folder.name="Movement" )
      movement = ss [ which( is.finite(rowSums( ss[,c(coords0, coords1)]))) , ]

      for ( i in 1:nrow(movement) ) { # low defined at top
        start = movement[i, coords0]
        end = movement[i, coords1] 
        names(end) = names(start)
        lc = rbind( start, end )
        kml.line.make( con, style.id="line.movement", x=lc  )
      }
    kml.folder.end (con)
        
  kml.end( con )
  print( outfile )




