
make.mortality = function(p, redo=F) {
  # est total mortality
  outfilename = file.path( p$annual.results, "mortality.rates.rdata" )
  if (!redo) {
    load( outfilename )
    return ( xx )
  }

  save( xx, file=outfilename, compress=T )
  return(xx)
} 


