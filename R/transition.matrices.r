
  transition.matrices = function(p, sex=male, redo=F) {
    # a wrapper function to obtain transition matrices
    outfilename = file.path( p$annual.results, "transition.matrices.rdata" )
    if (!redo) {
      load( outfilename )
      return ( tmatrix )
    }

    # initialize matrices
    eps = 1e-4  # assume this is min resolution
    xTM =  array(data=NA, dim=c(p$nnodes,p$nnodes,p$nyears,p$nregions), dimnames=list(p$nodes,p$nodes,p$nodeyears,p$regions) )
    xFM =  array(data=NA, dim=c(p$nnodes,p$nyears,p$nregions), dimnames=list(p$nodes,p$nodeyears,p$regions) )
    
    # create transition matrices
    for (region in p$regions) {
      # observer based landings broken down by size-mat classes
      fm = landings.decomposed ( p, sex=sex, region )
      # kriged numerical abundance estimates
      ts = get.annual.timeseries( p, sex=sex, outtype="yearclass", region=region)
      k = ts$TS
      k[ which(k < eps) ] = 0
      for (y in p$nodeyears) {

        nnodes = length(p$nodes)
        yc0 = as.character(y-1)
        yc1 = as.character(y)

        # the transition matrix estimates the transfer function between
        # the post-fishery numbers [N.post(t-1)] to the prefishery Numbers at time t [N.pre(t)]
        # == N.pre[t] / N.post[t-1]
        # k contains N.post from kriged estimates
        # fm contains N.fishing from observer data
        # N.pre[t] = N.post[t] + N.fishing[t]
        # N.post[t] must be offset to be N.post[t-1] for transition matrix calculations

        # fall surveys from 2003 to present 
        # spring surveys from 1998 to 2001
        # spring and fall surveys in 2002
        # correct timing of landings to account for change in survey season (by offsetting landings)
        fm[,,"2002"] = fm[,,"2002"] + fm[,,"2001"]
        for (y in 2001:1999)  fm[,,as.character(y)] = fm[,,as.character(y-1)]
        fm[,,"1998"] = fm[,,"1998"] * NA

        N.pre  = k + fm  
        N.post1 = offset.ageclass(k)
        N.post0 = k  # no offsets

        TM = array( data=0, c(nnodes, nnodes), dimnames=list(p$nodes, p$nodes) )
        FM = array( data=0, c(nnodes), dimnames=list(p$nodes) )

        if (sex==male) {

          # imm->imm
          ii = N.pre[,"imm", yc1] / N.post1[,"imm",yc0]  
          ii = ifelse( !is.finite(ii) | ii < eps,  mean( ii[ ii > eps &  is.finite(ii) ]) , ii)

          # imm->sm; no offset in size as they do not change in size ("instar")
          sm = N.pre[,"imm.sm",yc1] / N.post0[,"imm",yc0]  
          sm = ifelse( sm < eps | !is.finite(sm), mean( sm[ sm > eps & is.finite(sm) ]), sm)

          # fraction of imm -> imm.sm
          fraction.sm.imm = N.post1[,"imm.sm",yc0] / ( N.post1[,"imm",yc0] + N.post1[,"imm.sm",yc0] )
          fraction.sm.imm = ifelse( fraction.sm.imm <= eps |  fraction.sm.imm >= 1 | !is.finite( fraction.sm.imm  ), 
            mean( fraction.sm.imm [ which( fraction.sm.imm > eps & fraction.sm.imm < 1 & is.finite( fraction.sm.imm ) ) ]), 
            fraction.sm.imm  )

          # imm->cc12  (assuming some fraction came from imm)
          is = (1 - fraction.sm.imm) * N.pre[,"CC1to2",yc1] / N.post1[,"imm",yc0]  
          is = ifelse( (is <= eps |  is >= 1 | !is.finite(is) ), mean( is[ is > eps & is < 1 & is.finite(is) ]), is)

          # sm -> CC12 (" ..."  sm)
          m0 = fraction.sm.imm * N.pre[,"CC1to2",yc1] / N.post1 [,"imm.sm",yc0]  
          m0 = ifelse( (m0 <= eps | !is.finite(m0) ),  mean( m0[ m0 > eps & is.finite(m0) ]), m0)
           
          # assume mortality 0.33 = 1/3yr (expected duration of the CC2to4 stage)
          expected.longevity.CC34.yr = 3 # assuming pseudo-steady-state (ie. 1/3 = 33% annual mortality)
          
          # the mortality percent per year is expressed as passage into the CC5 category
          loss.rate.CC34 = 1 / expected.longevity.CC34.yr  
          transfer.rate.CC34 = 1 - loss.rate.CC34
          m2 = rep(transfer.rate.CC34, length(ii))  # cc34 -> cc34 ; assume unit catchability
          names(m2) = names(ii)

          # this is the fraction that is assumed to have come from CC12
          difference = N.pre[,"CC3to4",yc1] - (transfer.rate.CC34 * N.post0[,"CC3to4",yc0])  
          phi = difference / N.post0[,"CC1to2",yc0] # cc12 -> cc34
          phi = ifelse( phi <= eps |  !is.finite(phi), mean( phi[ phi > eps & is.finite(phi) ]), phi)

          # cc34 -> cc5 
          mt = rep(loss.rate.CC34, length(ii) )
          names(mt) = names(ii)

          TM["imm.6",   "imm.5"] = ifelse( (ii["6"]==0 | !is.finite(ii["6"])), 1, ii["6"] ) # no fishing mort and undersampled .. assume 1 to be safe if est is too low
          TM["imm.7",   "imm.6"] = ifelse( (ii["7"] ==0 | !is.finite(ii["7"])), 1, ii["7"] )
          TM["imm.8",   "imm.7"] = ifelse( (ii["8"] ==0 | !is.finite(ii["8"])), 1, ii["8"] )
          TM["imm.9",   "imm.8"] = ii["9"]
          TM["imm.10",  "imm.9"] = ii["10"]
          TM["imm.11", "imm.10"] = ii["11"]
          TM["imm.12", "imm.11"] = ii["12"]

          TM["CC1to2.9",   "imm.8"] = is["9"]
          TM["CC1to2.10",  "imm.9"] = is["10"]
          TM["CC1to2.11", "imm.10"] = is["11"]
          TM["CC1to2.12", "imm.11"] = is["12"]

          TM["imm.sm.9",   "imm.9"] = sm["9"]
          TM["imm.sm.10",  "imm.10"] = sm["10"]
          TM["imm.sm.11",  "imm.11"] = sm["11"]
          TM["imm.sm.12",  "imm.12"] = sm["12"]

          TM["CC1to2.10", "imm.sm.9"] = m0["10"]
          TM["CC1to2.11","imm.sm.10"] = m0["11"]
          TM["CC1to2.12","imm.sm.11"] = m0["12"]
          TM["CC1to2.13","imm.sm.12"] = m0["13"]

          TM["CC3to4.9", "CC1to2.9"] = phi["9"]
          TM["CC3to4.10", "CC1to2.10"] = phi["10"]
          TM["CC3to4.11", "CC1to2.11"] = phi["11"]
          TM["CC3to4.12", "CC1to2.12"] = phi["12"]
          TM["CC3to4.13", "CC1to2.13"] = phi["13"]

          TM["CC3to4.9", "CC3to4.9"] = m2["9"]
          TM["CC3to4.10", "CC3to4.10"] =  m2["10"]
          TM["CC3to4.11", "CC3to4.11"] =  m2["11"]
          TM["CC3to4.12", "CC3to4.12"] =  m2["12"]
          TM["CC3to4.13", "CC3to4.13"] =  m2["13"]

          TM["CC5.9", "CC3to4.9"] = mt["9"]
          TM["CC5.10", "CC3to4.10"] = mt["10"]
          TM["CC5.11", "CC3to4.11"] = mt["11"]
          TM["CC5.12", "CC3to4.12"] = mt["12"]
          TM["CC5.13", "CC3to4.13"] = mt["13"]

        # Fishing mortality

          f.ii = fm[,"imm", yc1] / N.pre[,"imm",yc1]  # exploitation of imm
          f.ii = ifelse(  f.ii <= eps | f.ii >=1  | !is.finite(f.ii), mean( f.ii[ f.ii > eps & f.ii < 1 & is.finite(f.ii) ]), f.ii)
          
          f.sm = fm[,"imm.sm", yc1] / N.pre[,"imm.sm",yc1]
          f.sm = ifelse(  f.sm <= eps | f.sm >=1  | !is.finite(f.sm), mean( f.sm[ f.sm > eps & f.sm < 1 & is.finite(f.sm) ]), f.sm)
          
          f.is = fm[,"CC1to2", yc1] / N.pre[,"CC1to2",yc1]
          f.is = ifelse(  f.is <= eps  | f.is >=1 | !is.finite(f.is), mean( f.is[ f.is > eps & f.is < 1 & is.finite(f.is) ]), f.is)

          f.phi = fm[,"CC3to4", yc1] / N.pre[,"CC3to4",yc1]
          f.phi = ifelse(  f.phi <= eps  | f.phi >=1 | !is.finite(f.phi), mean( f.phi[ f.phi > eps & f.phi < 1 & is.finite(f.phi) ]), f.phi)

          f.mt = fm[,"CC5", yc1] / N.pre[,"CC5",yc1]
          f.mt = ifelse(  f.mt <= eps  | f.mt >=1 | !is.finite(f.mt), mean( f.mt[ f.mt > eps & f.mt < 1 & is.finite(f.mt) ]), f.mt)

          FM["imm.10"] = f.ii["10"]
          FM["imm.11"] = f.ii["11"]
          FM["imm.12"] = f.ii["12"]

          FM["imm.sm.10"] = f.sm["10"]
          FM["imm.sm.11"] = f.sm["11"]
          FM["imm.sm.12"] = f.sm["12"]

          FM["CC1to2.10"] = f.is["10"]
          FM["CC1to2.11"] = f.is["11"]
          FM["CC1to2.12"] = f.is["12"]
          FM["CC1to2.13"] = f.is["13"]

          FM["CC3to4.10"] = f.phi["10"]
          FM["CC3to4.11"] = f.phi["11"]
          FM["CC3to4.12"] = f.phi["12"]
          FM["CC3to4.13"] = f.phi["13"]

          FM["CC5.10"] = f.mt["10"]
          FM["CC5.11"] = f.mt["11"]
          FM["CC5.12"] = f.mt["12"]
          FM["CC5.13"] = f.mt["13"]
        
        }
            
        xTM[,,as.character(y),region] = TM
        xFM[,as.character(y),region] = FM
      }
    } #end regions

    tmatrix = list(TM=xTM, FM=xFM)
    save( tmatrix, file=outfilename, compress=T )
    return(tmatrix)
  } 



