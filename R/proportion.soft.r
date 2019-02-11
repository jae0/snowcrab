
proportion.soft = function(odb, region, year) {

  # sex codes
    male = 0
    female = 1
    sex.unknown = 2


  out= NULL

# Remove CW's outside norms and remove production (pre-sorted) samples
  i = which( odb$sex==male & odb$fishyr==year )  # --- fishing year is used and not the actual year caught
  r = aegis::polygon_inside(x=odb, region=aegis::polygon_internal_code(region), planar=FALSE)

  c.r.y = unique( intersect (r, i) )

  soft = which( odb$durometer < 68 )

  nsoft = length( intersect( c.r.y, soft ) )

  ntot = length( c.r.y )

  out = c( round(nsoft/ntot*100,2), nsoft, ntot)

  return(out)
}



