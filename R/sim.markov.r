
  sim.markov = function( ip=NULL, p, DS="file" ) {
      
    if (exists( "init.files", p)) LoadFiles( p$init.files ) 
    if (exists( "libs", p)) RLibrary( p$libs ) 
    if (is.null(ip)) ip = 1:p$nruns
   
    fn.root =  file.path( project.datadirectory("snowcrab"), "markov"  )
    fn =  file.path( project.datadirectory("snowcrab"), "markov", paste( p$senario, "rdata", sep=".")  )

    dir.create( fn.root, recursive=T, showWarnings=F )

    if ( DS=="saved" ) { # default is to return the saved output
      load (fn )
      return(fp)
    }

    if ( DS=="collect" ) {
      fb = proj.summ = NULL
      fl = file.path( fn.root, paste( "markov.sim", p$scenario,  ip, "rdata", sep=".") )
      for ( i in fl ) {
        load (fl )
        fb = rbind( fb, pp.fp$FB )
        proj.summ = rbind( proj.summ, pp.fp$proj.summ)
      } 
      fp = list( FB=fb, proj.summ=proj.summ )
      save( fp, file=fn, compress=T )
      return ( fp )
    }

 
    prj = FB = out =NULL
    
    for ( i in ip ) {
      ER = p$runs[, "er"]
      sizes = mean.weights.by.category( p )  # use the previously saved file
      
      tmatrix = transition.matrices(p) 
      tmatrix$TM = tmatrix$TM[,, p$good.years ,] # good TM estimates
       
      CY = as.character(p$start.projection.year)

      years.projection = c( p$start.projection.year: (p$start.projection.year + p$nyears.projection))
      
      eps = 1e-6
      N.dimnames = list( p$nodes, as.character(c( p$start.projection.year : (p$start.projection.year+p$nyears.projection))), p$regions )
      FB.dimnames = list( as.character(c(p$start.projection.year : (p$start.projection.year+p$nyears.projection))), 
                          c("fb.prefishery", "error", "landings","rate", "fb.postfishery", "fb.increment"), p$regions )
      
      N0 = array( data=NA, dim=c(p$nnodes, p$nyears.projection+1, p$nregions), dimnames=N.dimnames )
      ERp = B0 = B0sd = B1 = B1sd = L = Lsd = N1 = N1sd = N0sd = N0 
      FB = array( data=NA, dim=c(p$nyears.projection+1, 6, p$nregions), dimnames=FB.dimnames )
         
      for (region in p$regions) {
         
        # kriged numerical abundance estimates
          ts = get.annual.timeseries( p, sex=male, outtype="yearclass", region=region)
          # "ts$sd" is actually 2*SD .. the 95% CI
          Nk = ts$TS
          Nk[ Nk < eps] = 0
          Nk = Nk[ ,, as.character(p$nodeyears)]
      
          fm = tmatrix$FM
          fm[!is.finite(fm)] = 0
          
          tm = choose.transition.matrix( tmatrix$TM, region, threshold.to.delete=2, use.global.average=T ) # threshold.to.delete is the upper bound of a TM to believe 
          XM = tm$XM
          XMsd = tm$XMsd
          tm = NULL

          fishing = landings.decomposed ( p, sex=male, region )
          # fall surveys from 2003 to present 
          # spring surveys from 1998 to 2001
          # spring and fall surveys in 2002
          # correct timing of landings to account for change in survey season (by offsetting landings)
          fishing[,,"2002"] = fishing[,,"2002"] + fishing[,,"2001"]
          fishing[,,"1998"] = fishing[,,"1998"] * NA
          for (y in c(2001:1999)) {
            fishing[,,as.character(y)] = fishing[,,as.character(y-1)]
          }
          fishing = fishing[,,as.character(p$nodeyears)]

          Nall = Nk + fishing

          for (iy in 1:(p$nyears.projection+1)) {
            
            b = Lb = NULL
            
            if (iy==1) { 

              # prefishery data .. seed projection matrices with p$start.projection.year results
              N0[,iy,region] = convert.to.vector( Nall[,,CY], p$nodes) # number estimate (pre-fishery), (in kn)
              N0sd[,iy,region] = convert.to.vector( ts$TS.sd[,,CY], p$nodes) # ts$TS.sd is already 2 SD in number 

              b = estimate.biomass( N0[,iy,region], N0sd[,iy,region], sizes, y=p$start.projection.year, r=region)  
              B0[,iy,region] = b$x
              B0sd[,iy,region] = b$error
              
              FB[iy,c("fb.prefishery","error"),region] = subset.biomass (B0[,iy,region], B0sd[,iy,region], p$nodes, type=p$scenario)  

              fishery.stats = get.fishery.stats.by.region(Reg=region, y=p$fisheryyears )
              ERp[,iy,region] = fm[, CY, region]

              # landings expected based upon exploitation strategy 
              L[,iy,region] = N0[,iy,region] * ERp[,iy,region]
              Lsd[,iy,region] = L[,iy,region] * N0sd[,iy,region]/N0[,iy,region]
              Lb = estimate.biomass (L[,iy,region], Lsd[,iy,region], sizes, y=p$start.projection.year, r=region) # use last known distribution of sizes
              FB[iy,"landings",region] = subset.biomass (Lb$x, Lb$error, p$nodes, type=p$scenario)[1]

              # update N and B vectors
              B1[,iy,region]   = B0[,iy,region] - Lb$x
              B1sd[,iy,region] = error.propagate( error.x=B0sd[,iy,region], error.y=Lb$error, type="sum" )
          
              N1[,iy,region]   = N0[,iy,region] - L[,iy,region]  
              N1sd[,iy,region] = error.propagate( error.x=N0sd[,iy,region], error.y=Lsd[,iy,region], type="sum" )
              
              FB[iy,"fb.increment", region] =  NA
              FB[iy,"fb.postfishery",region] = FB[iy,"fb.prefishery",region] - FB[iy,"landings",region] 
              FB[iy,"rate",region] = FB[iy,"landings",region] / (FB[iy,"fb.prefishery",region] ) 
             
            } else {
              
              # forward project the prefishery numbers
              N0[,iy,region] = XM %*% N1[,(iy-1),region] # projection (pre-fishery)
              N0sd[,iy,region] = N1[,(iy-1), region] * error.propagate(x=N1[,(iy-1),region], error.x=N1sd[,(iy-1),region], y=XM, error.y=XMsd, type="product"  )

              # projected fishable biomass estimates broken down by class (pre-fishery)
              b = estimate.biomass( N0[,iy,region], N0sd[,iy,region], sizes, y=p$start.projection.year, r=region)  # use  p$start.projection.year size data where possible
              B0[,iy,region] = b$x
              B0sd[,iy,region] = b$error

              FB[iy,c("fb.prefishery","error"),region] = subset.biomass (B0[,iy,region], B0sd[,iy,region], p$nodes, type=p$scenario) 

              target.landings = FB[iy,"fb.prefishery",region] * ER
              ER.projected = projected.exploitation.rate( target.landings, fm.reg=fm[,,region], 
                  numbers.initial=N0[,iy,region], sizes=sizes, yr=p$start.projection.year, region=region, scenario=p$scenario )
              ERp[,iy,region] = ER.projected$exploitation.rate.optim
              
              # landings expected based upon exploitation strategy 
              L[,iy,region] = ER.projected$numbers.optim 
              Lsd[,iy,region] = L[,iy,region] * N0sd[,iy,region]/N0[,iy,region]  # linear scaling of error
              Lb = estimate.biomass (L[,iy,region], Lsd[,iy,region], sizes, y=p$start.projection.year, r=region) # use last known distribution of sizes
              FB[iy,"landings",region] = subset.biomass (Lb$x, Lb$error, p$nodes, type=p$scenario)[1]
            
            # update N and B vectors
              B1[,iy,region]   = B0[,iy,region] - Lb$x
              B1sd[,iy,region] = error.propagate( error.x=B0sd[,iy,region], error.y=Lb$error, type="sum" )
          
              N1[,iy,region]   = N0[,iy,region] - L[,iy,region]  
              N1sd[,iy,region] = error.propagate( error.x=N0sd[,iy,region], error.y=Lsd[,iy,region], type="sum" )
              
              FB[iy,"fb.increment", region] =  FB[iy,"fb.prefishery",region] - (FB[(iy-1),"fb.postfishery",region] )
              FB[iy,"fb.postfishery",region] = FB[iy,"fb.prefishery",region] - FB[iy,"landings",region] 
              FB[iy,"rate",region] = FB[iy,"landings",region] / (FB[iy,"fb.prefishery",region] ) 
              
            }
            
          } # end for projection year
         

        } # end for region
        
        FB = round(FB, digits=5)
     
        fb.projection = list(FB=FB, N0=N0, N0sd=N0sd, B0=B0, B0sd=B0sd, N1=N1, N1sd=N1sd, B1=B1, B1sd=B1sd, ERp=ERp, L=L ) 
     
        FB = rbind( FB, cbind(as.data.frame.table(fb.projection$FB), ER=ER ) )
        names(FB) = c("yr", "vars", "region", "val", "ER")
        FB = factor2number(FB, c("yr", "val", "ER") )

        out = rbind( out, cbind(as.data.frame.table(fb.projection$N0), ER=ER, datatype="N0" ) )
        out = rbind( out, cbind(as.data.frame.table(fb.projection$N0sd), ER=ER, datatype="N0sd" ) )
        out = rbind( out, cbind(as.data.frame.table(fb.projection$N1), ER=ER, datatype="N1" ) )
        out = rbind( out, cbind(as.data.frame.table(fb.projection$N1sd), ER=ER, datatype="N1sd" ) )
        out = rbind( out, cbind(as.data.frame.table(fb.projection$B0), ER=ER, datatype="B0" ) )
        out = rbind( out, cbind(as.data.frame.table(fb.projection$B0sd), ER=ER, datatype="B0sd" ) )
        out = rbind( out, cbind(as.data.frame.table(fb.projection$B1), ER=ER, datatype="B1" ) )
        out = rbind( out, cbind(as.data.frame.table(fb.projection$B1sd), ER=ER, datatype="B1sd" ) )
        out = rbind( out, cbind(as.data.frame.table(fb.projection$ERp), ER=ER, datatype="ERp" ) )
        out = rbind( out, cbind(as.data.frame.table(fb.projection$L), ER=ER, datatype="L" ) )
          

        names(out) = c("node", "yr", "region", "val", "ER", "datatype")
        proj.summ = factor2number(out, c("yr", "val", "ER") )

        prj = list( FB=FB, proj.summ=proj.summ )
        
        save ( prj, file = file.path(  fn.root, paste( "markov.sim", p$scenario, i, "rdata", sep="." )) ,compress=T ) 
        print ( ER )
      }
      return ("Completed")
  }

    

