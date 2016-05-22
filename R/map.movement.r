
 map.movement = function( p, outdir ) {
 
    tags.datadir= file.path( project.datadirectory("snowcrab"), "data", "tagging" )

    # plot movement of crab in maps
    # and generated interpolated fields of movemnent magnitudes, orientations and their errors

    marked.file = "tags.1996_2001.csv"
    marked2.file = "tags_summary1993_2005.csv"
    recaps.file = "recaptures.csv"


    move = get.move (redo=T, region="cfaall") #  

    toplot = move
    toplot = toplot[ is.finite(toplot$lon0 + toplot$lat0 + toplot$lon1 + toplot$lat1) ,]
    toplot = toplot[, c("lon0", "lat0", "lon1", "lat1")]

    p$arrow = "-Svs0.01c/0.05c/0.03c -G50/50/250"

    p$basedir = outdir
    p$outfile.basename = file.path(p$basedir, "arrows" )

    if (redo.basemap) gmt.basemap (p)

    gmt.arrows(p, toplot, conversions=p$conversions )
   
    pause(30)
    files = sort( list.files( p$basedir, filter="[*.ps]$", all.files=T, full.names=T, recursive=T ))
    remove.files ( files ) 

  }

