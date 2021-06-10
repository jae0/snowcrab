map.fisheries.data = function(
  outdir, 
  yrs, 
  variables=c('effort', 'cpue', 'landings') ,
  FUN = c(sum, mean, sum),
  probs=c(0.025, 0.975), 
  pres = 10,
  plot_crs = st_crs( projection_proj4string("lonlat_wgs84"))  ,
  outformat="pdf",
  ...
  ) {

  require(sf)
  require(tmap)
  
  x = logbook.db( DS="logbook" )
  x = x [polygon_inside( x, region="isobath1000m"),]
  x = x[ which(x$effort <= 300) ,]
  x = x[ which(x$cpue < 500),]
  x$year = x$yr #this creates proper offset for 4X, 2017-18 season =2017
  x$landings = x$landings/1000 

  
  x = st_as_sf( x, coords= c("lon", "lat"), crs=st_crs( projection_proj4string("lonlat_wgs84") ) )
  x = st_transform(x, st_crs(p$aegis_proj4string_planar_km ) )
  
  if (missing(yrs)) yrs=sort(unique(x$year))

  sppoly = st_as_sf( st_make_grid( x, cellsize=pres, what="polygons", square=TRUE ) )
  sppoly$AUID = as.character( 1:nrow(sppoly) )  # row index
  row.names(sppoly) = sppoly$AUID

  x$AUID = st_points_in_polygons( pts=x, polys=sppoly[, "AUID"], varname="AUID" )


  coastline = aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs ) 
  
  # depth contours
  isobaths = aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs  )

  managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
  
  ellps = list(...)

  for (v in 1: length(variables)) {
    vn =variables[v]
    er = quantile( x[[vn]], probs=probs, na.rm=TRUE )
    datarange = pretty( er  )  
    
    outloc = file.path( outdir, vn )
    dir.create (outloc, showWarnings=FALSE, recursive =TRUE)

    for (i in 1:length(yrs)){

      outfn = paste( vn, yrs[i], sep=".")
      iy = which(x$year == yrs[i] )

      toplot = x[[vn]] [iy]  

      oo = tapply( toplot, x[["AUID"]] [iy], FUN[[v]], na.rm=TRUE )
 
      sppoly$z = NA
      sppoly$z[ match( names(oo) , sppoly$AUID ) ] = oo

      sppoly$z[ sppoly$z < er[1] ] = er[1] # set levels higher than max datarange to max datarange
      sppoly$z[ sppoly$z > er[2] ] = er[2] # set levels higher than max datarange to max datarange

      # do the map
      fn = file.path( outloc, paste(outfn, outformat, sep="." ) )

      if (outformat=="pdf") pdf( file=fn, width=8, height=6, bg='white', pointsize=20 )
      if (outformat=="svg") svg( filename=fn, width=8, height=6, bg='white', pointsize=20   )
      if (outformat=="png") png( filename=fn, width=1200, height=800, pointsize=20, res=100 )

        tmap_mode("plot")
        o = tm_shape( sppoly, projection=plot_crs ) +
          tm_polygons(
            "z",
            style = ifelse ( exists("style", ellps), ellps[["style"]], "cont" ) ,
            breaks = datarange,
            title= paste( vn, ":", yrs[i] ),
            border.col = NULL,
            colorNA = NULL,
#            constrast=c(0,0.6),
            showNA=FALSE,
            lwd = 0.5, 
            palette = ifelse ( exists("palette", ellps), ellps[["palette"]], "YlOrRd"),
            border.alpha = 0.5,
            alpha = 0.95,
            legend.is.portrait = FALSE ) +
        tm_shape( coastline, projection=plot_crs ) +
          tm_polygons( col="grey80" ) +
        tm_shape( isobaths, projection=plot_crs ) +
          tm_lines( col="lightgray", alpha=0.5) +
        tm_shape( managementlines, projection=plot_crs ) +
          tm_lines( col="grey40", alpha=0.6, lwd=2) +

        tm_compass( position=c( "right", "top")) + 
        tm_scale_bar( position=c("right", "bottom" ), width=0.2, text.size=0.7) +
        tm_legend( position=c("left", "top") ,  frame=TRUE, scale = 1 , title.size=1.5, text.size=0.80, legend.width=0.75) +
        tm_layout( frame=FALSE )

        print(o)
      dev.off()
      print(fn)
    }

  }

}
