
  map.cat.information = function(p, outdir ) {
    bcat = snowcrab.db( DS="cat.georeferenced")
    species = c( "cod", "skates", "thornyskate", "northernshrimp",
               "brittlestar", "basketstar",
               "snowcrab",  "hermitcrab", "jonahcrab",
               "lessertoadcrab", "northernstonecrab", "porcupinestonecrab", "portlyspidercrab",
               "redcrab", "atlanticrockcrab", "toadcrab", "snowcrab",'atlanticwolffish','halibut','smoothskate','winterskate')
    p$tension = "-T.4"  # 0.35+ for steep; 0.25 for smooth
    p$maskres = "-S16k"
    p$interpres = "-nb"
    
    
    for (sp in species) {
      outvars = c("trip", "set", "yr", "lon", "lat", "totmass", "totno", "sa")
      basedir = file.path( outdir, p$spatial.domain, sp )
      sset = bcat[ taxonomy.filter.taxa (bcat$spec, taxafilter=sp, outtype="groundfishcodes" ), outvars]
      if (dim(sset)[1] > 1) {
        variables = c("totno")
        sset = sset[ is.finite(rowSums(sset[, c("yr", "lon", "lat")])), ]
        sset$totno = log10( sset$totno )  ###### ---------------------- log transform done here as "totno" is not part of the snowcrab recoed db
        # must  sum the catches or numbers by set/trip (across species) before plotting
        # incase multiple species are selected
        gmt.map.variables( sset, params=p, variables=variables, plottimes=p$plottimes, basedir=basedir, conversions=p$conversions, init.files=p$init.files )
       }
    }
    return( "Done" )
  }


