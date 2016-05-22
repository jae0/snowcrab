
  fishing.area.designations = function( x, type="lonlat" ) {

    if (type=="lonlat") planar=F
    if (type=="planar") planar=T

    icfa4x = filter.region.polygon( x, "cfa4x", planar=planar)
    icfanorth = filter.region.polygon( x, "cfanorth", planar=planar)
    icfasouth = filter.region.polygon( x, "cfasouth", planar=planar)
    G = rep( NA, nrow( x ) )
    G[icfa4x] = "cfa4x"
    G[icfanorth] = "cfanorth"
    G[icfasouth] = "cfasouth"
    x$cfa = G
    x$cfa.factor = factor( x$cfa, levels=c("cfa4x","cfanorth" ,"cfasouth") )
    return(x)
  }
    

