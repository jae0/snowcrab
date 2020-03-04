p = bio.snowcrab::snowcrab_carstm( DS="parameters",
                                   assessment.years=1999:2019 )

bio = snowcrab_carstm(p=p, DS="carstm_output_spacetime_biomass",carstm_model_label=p$carstm_model_label  )
num = snowcrab_carstm(p=p, DS="carstm_output_spacetime_number", carstm_model_label=p$carstm_model_label  )

u = which.max( bio[,"2019"] )
dimnames(bio)
# givem the AUID and years
dimnames(bio)$AUID[ u]

M = snowcrab_carstm( p=p, DS="carstm_inputs" )  # will redo if not found M[M$AUID==dimnames(bio)$AUID[ u] & M$year==2019 ,]



# to plot:

sppoly = areal_units( p=p )  # to reload

vn = "pred"
sppoly@data[,vn] = bio[,"2019"] # adds 2019 data to polygons

brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette),
                      style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )

sppoly@data[u,vn] = 10000  # force to go out of range

spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


#determine which station ID's from survey are associated with which AU's for the model run?

#see line 215 of snowcrab_carstm.R

crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))
sppoly = areal_units( p=p )  # will redo if not found

M$AUID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area







#--------------------------------------------------------------------------------

N = snowcrab.db( p=p, DS="biological_data" )

 N$density = N$totno /  N$data_offset

hist( N$density  )

#compare to:

hist( N$density[ which(N$yr==2019)] )
> > > > > >
  > > >
  > > >

