
  landings.decomposed = function( p, sex, REGION) {
    
    # observer data
    odb = observer.db("odb")

    sizes = mean.weights.by.category( p )  # use the previously saved file

    catch.odb = observer.get.counts.by.class( p, odb[filter.region.polygon(odb,recode.areas( REGION )) ,] )
    cl = make.classes(sex)

    # obtain landings
    fishery.stats = get.fishery.stats.by.region( Reg=REGION, y=p$fisheryyears )

    fishing = array( data=0, dim=c(length(cl$yclass), length(cl$cats), length(p$fisheryyears)),
                    dimnames=list(cl$yclass, cl$cats, p$fisheryyears) )

    for (jj in p$fisheryyears) {
      kk = which(catch.odb$yr==jj)
      if ( length(kk) == 0 ) next()
      catch.odb.yr = catch.odb[kk, ]
      catch.N = convert.to.vector(catch.odb.yr, p$nodes, type="observer.data" )
      catch.bio.est = estimate.biomass (catch.N, catch.N*0.01, sizes, y=jj, r=REGION)
      catch.B = catch.bio.est$x
      total.observed.biomass =  sum(catch.B, na.rm=T) 

      catch.rescaled =  catch.B /total.observed.biomass
    
      catch.total =  catch.rescaled * fishery.stats$landings[ which(fishery.stats$yr==jj) ] * 1000 # 1000 is to convert kg to tons
      FM.N = estimate.number( catch.total, catch.bio.est$error, sizes, y=jj, r=REGION)
      fishing[,,as.character(jj)] = nodes2matrix (FM.N$x, cl)
    }

    return (fishing)
  } 


