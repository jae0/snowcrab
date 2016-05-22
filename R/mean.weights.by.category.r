
  mean.weights.by.category = function(p, redo=F) {

    outfilename = file.path( p$annual.results, "crab.sizes.rdata" )

    if (!redo) {
      load( outfilename )
      return( sizes )
    }

    set = snowcrab.db("set.merge.det")
    set = set[,c("trip", "set", "yr", "lon", "lat", "sa")]
     
    det = snowcrab.db("det.initial")
    det = det[,c("trip", "set", "sex", "cw", "mass", "abdomen", "chela", "mat", "shell", "gonad", "eggcol", "eggPr", "durometer" )]
    det = merge(det, set, by=c("trip","set"), all.x=T, all.y=F, sort=F)
    rm(set)

    p$regions = c(p$regions, "cfaall" ) # add a global category as well
    p$nregions = length(p$regions)
    
    results = c("mass.mean", "mass.sd", "cw.mean", "cw.sd", "n")
    nresults = length(results)
    
    sizes = array(data=NA, dim=c(p$nnodes,nresults,p$nfisheryyears,p$nregions), 
      dimnames=list(p$nodes, results, p$fisheryyears, p$regions)  )

    det.years = det$yr
    det.locs = det[,c("lon","lat")]

    for (region in p$regions) {
    for (y in p$fisheryyears) {
    for (iv in 1:length(p$nodes)) {
      

      node = p$nodes[iv]
      nodelabel = p$nodelabels[iv] 
      yc = as.character(y)
      ic = filter.class( det, type=nodelabel )
      ir = filter.region.polygon(det.locs, region=region)
      iy = which( det.years==y )
      i = sort( intersect( intersect(  ir, iy) , ic )  )
      if (length(i)<1) next
      sizes[ node, "mass.mean", yc, region ] = weighted.mean(x=det$mass[i], w=det$sa[i], na.rm=T)
      sizes[ node, "mass.sd", yc, region ] = sd(x=det$mass[i], na.rm=T)
      sizes[ node, "cw.mean", yc, region ] = weighted.mean(x=det$cw[i], w=det$sa[i], na.rm=T)
      sizes[ node, "cw.sd", yc, region ] = sd(x=det$cw[i], na.rm=T)
      sizes[ node, "n", yc, region ] =  length(i)
    
    }}}
   
    # fill in from each area 
    for (region in p$regions) {
    for (iv in 1:length(p$nodes)) {
      node = p$nodes[iv]
      for (v in results ) { 
        i = which(!is.finite( sizes[ node, v, , region ] ))
        sizes[ node, v, i,region ] = apply(sizes[ , v, , region ], 1, mean,na.rm=T) [iv]
      }
    }}
    
    # fill in using global means where possible
    for (region in p$regions) {
    for (iv in 1:length(p$nodes)) {
      node = p$nodes[iv]
      for (v in results ) { 
        i = which(!is.finite( sizes[ node, v, , region ] ))
        sizes[ node, v, i,region ] = apply(sizes[ , v, , 4 ], 1, mean,na.rm=T) [iv]
      }
    }}
    
    save (sizes, file=outfilename, compress=T)
    return (sizes)
  }


