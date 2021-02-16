map.fisheries.data = function(
  outdir, 
  yrs, 
  variables=c('effort', 'cpue', 'landings') ,
  FUN = c(sum, mean, sum),
  probs=c(0.025, 0.975), 
  pres = 10,
  crs = st_crs( projection_proj4string("lonlat_wgs84"))  ,
  outformat="pdf"
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


  coastline = aegis.coastline::coastline_db( p=p, DS="eastcoast_gadm", project_to=crs )
  
  # depth contours
  isobaths = aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=crs  )


  managementlines = aegis.polygons::area_lines.db( DS="cfa.regions" )
  managementlines = st_sfc( st_multilinestring( sapply( managementlines, as.matrix) ) )
  st_cast(managementlines, "MULTILINESTRING")
  st_crs(managementlines) = st_crs( projection_proj4string("lonlat_wgs84")) 
  


  for (v in 1: length(variables)) {
    vn =variables[v]
    er = quantile( x[[vn]], probs=probs, na.rm=TRUE )
    datarange = seq( er[1], er[2], length.out=7) 
    
    outloc = file.path( outdir, vn )
    dir.create (outloc, showWarnings=FALSE, recursive =TRUE)

    for (i in 1:length(yrs)){

      outfn = paste( vn, yrs[i], sep=".")
      iy = which(x$year == yrs[i] )

      toplot = x[[vn]] [iy]  
      auid = x[["AUID"]] [iy]
      oo = tapply( toplot, auid, FUN[[v]], na.rm=TRUE )
      oo[ oo< er[1] ] = er[1]
      oo[ oo> er[2] ] = er[2]
      sppoly$z = NA
      sppoly$z[ match( names(oo) , sppoly$AUID ) ] = oo

      sppoly$z[ sppoly$z > er[2] ] = er[2] # set levels higher than max datarange to max datarange

      # do the map
      fn = file.path( outloc, paste(outfn, outformat, sep="." ) )

      if (outformat=="pdf") pdf( file=fn, width=9, height=7, bg='white', pointsize=12 )
      if (outformat=="svg") svg( filename=fn, width=9, height=7, bg='white', pointsize=12   )
      if (outformat=="png") png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )

        tmap_mode("plot")
        o = tm_shape( sppoly, projection=crs ) +
          tm_polygons(
            "z",
            style = "cont",
            breaks = datarange,
            title= passte( vn, ":", yrs[i] ),
            border.col = NULL,
            colorNA = NULL,
#            constrast=c(0,0.6),
            showNA=FALSE,
            lwd = 0.5, 
            palette = "YlOrRd",
            border.alpha = 0.5,
            legend.is.portrait = FALSE ) +
        tm_shape( coastline, projection=crs ) +
          tm_polygons( col="grey80" ) +
        tm_shape( isobaths, projection=crs ) +
          tm_lines( col="lightgray", alpha=0.5) +
        tm_shape( managementlines, projection=crs ) +
          tm_lines( col="grey20", alpha=0.75, lwd=2) +

        tm_compass( position=c( "right", "top")) + 
        tm_scale_bar( position=c("left", "bottom" ) ) +
        tm_layout( frame=FALSE, legend.text.size= 0.7 )
        print(o)
      dev.off()
      print(fn)
    }

  }

}
