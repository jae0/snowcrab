
proportion.soft = function(odb, region, year) {

  out= NULL

# Remove CW's outside norms and remove production (pre-sorted) samples
  i = which( odb$sex==male & odb$fishyr==year )  # --- fishing year is used and not the actual year caught
  r = filter.region.polygon(x=odb, region=recode.areas(region), planar=F)

  c.r.y = unique( intersect (r, i) )

  soft = which( odb$durometer < 68 )

  nsoft = length( intersect( c.r.y, soft ) )

  ntot = length( c.r.y )

  out = c( round(nsoft/ntot*100,2), nsoft, ntot)

  return(out)
}



