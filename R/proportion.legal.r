
  proportion.legal = function(odb, region, year) {

    out= NULL

# Remove CW's outside norms and remove production (pre-sorted) samples
    i = which( odb$sex==male & odb$prodcd_id=="0" & odb$cw > 50 & odb$cw < 170 & odb$fishyr==year )
    j = which( odb$sex==male & odb$prodcd_id=="0" & odb$cw >= 95 & (odb$durometer >= 68 | odb$shell >=2) )

    r = filter.region.polygon(x=odb, region=recode.areas(region), planar=F)
    total.n = intersect (r, i)
    total.com = intersect(total.n, j)

    out = c( length(total.n), length(total.com), round( length(total.com)/length(total.n) * 100) )

    return(out)
  }


