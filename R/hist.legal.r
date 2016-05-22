
  hist.legal = function( odb, region, year, bks = seq( 90, 175, 5 )) { 
    out= NULL

    # Remove CW's outside norms and remove production (pre-sorted) samples
    i = which( odb$sex==male & odb$prodcd_id=="0" & odb$cw >= 95 & odb$cw < 170 & odb$fishyr==year & (odb$durometer >= 68 | odb$shell >=2 ) & odb$mat==1 )

    r = filter.region.polygon(x=odb, region=recode.areas(region), planar=F)
    z = intersect (r, i)
    
    hh = hist (odb$cw[z], breaks=bks, plot=F )
    
    out = hh$counts

  }
  

