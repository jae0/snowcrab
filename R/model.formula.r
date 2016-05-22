
    model.formula = function( V="default", spatial.knots=100 ) {
      out = NULL
      if ( V %in% c("default", "R0.mass") ) {
          # basic model -- less spatial density .. all covariates
           out = paste( ' Y ~ s(dyear, bs="cr" ) ', 
            ' + s(tmean) + s(dt.annual, bs="cr" ) + s(dt.seasonal, bs="cr" ) ', 
            ' + s(tamp, bs="cr" )+ s(wmin, bs="cr" ) ', 
            ' + s(z) + s(dZ, bs="cr" )  + s(substrate.mean, bs="cr" )', 
            ' + s(plon, plat, k=', spatial.knots, ', bs="tp", by=as.factor(yr) ) + as.factor(yr)'  )  }
        
         
      if ( V %in% c("default.2014") ) {
          # basic model -- less spatial density .. all covariates
           out = paste( ' Y ~ s(dyear, bs="ts" ) ', 
            ' + s(tmean) + s(dt.annual, bs="ts" ) + s(dt.seasonal, bs="ts" ) ', 
            ' + s(tamp, bs="ts" )+ s(wmin, bs="ts" ) ', 
            ' + s(z) + s(dZ, bs="ts" )  + s(substrate.mean, bs="ts" )     ', 
            ' + s(ca1, bs="ts" ) + s(ca2, bs="ts" ) ', 
            ' + s(Npred, bs="ts") + s(Z, bs="ts" ) ', 
            ' + s(smr, bs="ts" ) + s(A, bs="ts" ) + s(Ea, bs="ts" )    ', 
            ' + s(nss.shannon, bs="ts" ) + s(nss.rsquared, bs="ts" ) ', 
            ' + s(plon, plat, k=', spatial.knots, ', bs="tp", by=as.factor(yr) ) + as.factor(yr)'  ) 
      }
      
      if ( V %in% c("R0.mass.2014") ) {
        #  same as basic model but with higher spatial densities
          out = formula( Y ~ te(dyear, bs="cs" ) 
            + te(tmean,bs='cs') + te(dt.annual, bs="cs" ) + te(dt.seasonal, bs="cs" ) 
            + te(tamp, bs="cs" )+ te(wmin, bs="cs" ) 
            + te(z,bs='cs') + te(dZ, bs="cs" )  + te(substrate.mean, bs="cs" )             
            + te(ca1, bs="cs" ) + te(ca2, bs="cs" ) 
            + te(Npred, bs="cs") + te(Z, bs="cs" ) 
            + te(smr, bs="cs" ) + te(A, bs="cs" ) + te(Ea, bs="cs" )     
            + te(nss.shannon, bs="cs" ) + te(nss.rsquared, bs="cs" ) 
            + te(plon, plat, k=150, bs="cs", by=as.factor(yr) ) + as.factor(yr)  ) 
      }

			if ( V=="simple" ) {
        # simple model in case of failures in other attempts
        out = formula( Y ~  s(dyear, k=4, bs="ts" ) 
            + s(tmean) + s(z) + s(dZ, k=4, bs="ts" ) + s(substrate.mean, k=4, bs="ts" )             
            + s(plon, plat, bs="tp", by=as.factor(yr)) + as.factor(yr) ) 
      }

      if ( V=="totno.male.mat" ) {
        # simple model in case of failures in other attempts
        out = formula( Y ~  s(dyear, k=4, bs="ts" ) 
            + s(tmean) + s(z) + s(dZ, k=4, bs="ts" ) + s(substrate.mean, k=4, bs="ts" )             
            + s(plon, plat, bs="tp", by=as.factor(yr)) + as.factor(yr) ) 
      }

      if(V == 'simple.2015'){
         out = formula( Y ~  s(t, k = 5, bs = 'ds') +  
            + s(tmean.cl , k = 5, bs = 'ds') + s(z, k=4,bs = 'ds') + s(dZ, k=4, bs="ds" ) + s(substrate.mean, k=4, bs="ds" )             
            + s(plon, plat, k=100,bs="ds",by=as.factor(yr)) +   s(A,k=4,bs='ds')+  s(mass,k=4,bs='ds')+
            s(Z,k=4,bs='ds') + as.factor(yr) )# nonsig-- s(smr,k=4,bs='ds')+s(Ea,k=4,bs='ds')+ s(Npred,k=4,bs='ds')+s(qm,k=4,bs='ds')
      }

			if ( V=="R0.mass.2012" ) {
        # enough data for a high spatial resolution of data variability
        out = formula( Y ~ s(dyear, k=4, bs="ts" ) 
            + s(tmean) + s(dt.annual, k=4, bs="ts" ) + s(dt.seasonal, k=4, bs="ts" ) 
            + s(tamp, k=4, bs="ts" )+ s( wmin, k=4 , bs="ts" ) 
            + s(z) + s(dZ, k=4, bs="ts" )  + s(substrate.mean, k=4, bs="ts" )             
            + s(ca1, k=4, bs="ts" ) + s(ca2, k=4, bs="ts" ) 
            + s(Npred, k=4 , bs="ts") + s(Z, k=4, bs="ts" ) + s(smr, k=4, bs="ts" ) + s(mr, k=4, bs="ts" ) 
            + s(plon, plat, bs="ts", k=spatial.knots, by=as.factor(yr)) + as.factor(yr) ) 
      }
     
      if ( V=="R0.mass.old" ) {
        # enough data for a high spatial resolution of data variability
        out = formula( Y ~  s(dyear, k=4) + s( t-t.annual ) + s(tmean) + s(t.annual-tmean) + s( tamp)+ s( wmin) 
            + s( plon, plat, k=400, by=as.factor(yr) ) + s( z ) + s(substrate.mean ) + s(ddZ) + s(dZ) ) 
      }

      if ( V=="R0.mass.2011" ) {
        # enough data for a high spatial resolution of data variability
        out = formula( Y ~  s(yr) + s(dyear, k=4, bs="ts" ) 
            + s(tmean) + s(dt.annual, k=4, bs="ts" ) + s(dt.seasonal, k=4, bs="ts" ) 
            + s(tamp, k=4, bs="ts" )+ s( wmin, k=4 , bs="ts" ) 
            + s(z) + s(dZ, k=4, bs="ts" )  + s(substrate.mean, k=4, bs="ts" )             
            + s(ca1, k=4, bs="ts" ) + s(ca2, k=4, bs="ts" ) 
            + s(Npred, k=4 , bs="ts") + s(Z, k=4, bs="ts" ) + s(smr, k=4, bs="ts" ) + s(mr, k=4, bs="ts" ) 
            + s(plon, plat, k=250, bs="tp") ) 
      }
    
      if ( V=="R0.mass.environmentals.only" ) {
        # enough data for a high spatial resolution of data variability
        out = formula( Y ~ s(tmean) + s(dt.annual, k=4, bs="ts" ) + s(dt.seasonal, k=4, bs="ts" ) 
            + s(tamp, k=4, bs="ts" )+ s( wmin, k=4 , bs="ts" ) 
            + s(z) + s(dZ, k=4, bs="ts" ) + s(ddZ, k=4, bs="ts" )   + s(substrate.mean, k=4, bs="ts" )             
            + s(plon, plat, k=250, bs="tp", by=as.factor(yr)) + as.factor(yr) ) 
      }
    
      if ( V=="R0.mass.2009" ) {
        # enough data for a high spatial resolution of data variability
        out = formula( Y ~ s(yr, bs="ts" ) + s(dyear, bs="ts" ) 
            + s( dt.seasonal, bs="ts" ) + s(tmean, bs="ts" ) + s(dt.annual, bs="ts" ) 
            + s( tamp, bs="ts" ) + s( wmin, bs="ts"  ) 
            + s( z , bs="ts" ) + s(substrate.mean, bs="ts" ) + s(ddZ, bs="ts" ) + s(dZ, bs="ts"  ) 
            + s( plon, plat, bs="ts", k=200 )  ) 
      }
      print(out)
      out = as.formula( out )
      return ( out )
    }



