
  hist.legal = function( odb, region, year, bks = seq( 90, 175, 5 )) {
    out= NULL

    # sex codes
    male = 0
    female = 1
    sex.unknown = 2

    # maturity codes
    immature = 0
    mature = 1
    mat.unknown = 2


    # Remove CW's outside norms and remove production (pre-sorted) samples
    i = which( odb$sex==male & odb$prodcd_id=="0" & odb$cw >= 95 & odb$cw < 170 & odb$fishyr==year & (odb$durometer >= 68 | odb$shell >=2 ) & odb$mat==1 )

    r = polygon_inside(x=odb, region=aegis.polygons::polygon_internal_code(region), planar=F)
    z = intersect (r, i)

    hh = hist (odb$cw[z], breaks=bks, plot=F )

    out = hh$counts

  }


