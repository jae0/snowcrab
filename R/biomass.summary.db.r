
  biomass.summary.db = function( DS="complete", p=NULL ) {
    #browser()
    
    if (DS %in% c("complete", "complete.redo") ) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "bugs.rdata" ) 
      if ( DS == "complete") {
        out = NULL
        if (file.exists(fn)) load(fn)
        return( out )
      }
      
      L = biomass.summary.db( DS="L.redo", p=p  )  # must go first as part of biomass estimates
      B = biomass.summary.db( DS="B.redo", p=p )  # rename to avoid confusion below as B is also used
      B.sd = biomass.summary.db( DS="B.sd.redo", p=p )  # rename to avoid confusion below as B is also used
      R = biomass.summary.db( DS="R.redo", p=p  )  # rename to avoid confusion below as B is also used
      R.sd = biomass.summary.db( DS="R.sd.redo", p=p  )  # rename to avoid confusion below as B is also used
      
      # geometric means -- sd are on log scale
      Bg = biomass.summary.db( DS="B_geomean.redo", p=p  )
      Bg.sd = biomass.summary.db( DS="B_geomean.sd.redo", p=p  )
      Rg = biomass.summary.db( DS="R_geomean.redo", p=p  )
      Rg.sd = biomass.summary.db( DS="R_geomean.sd.redo", p=p  )
 

      # cfa4x have had no estimates prior to 2004
      spring = 1998:2003
      ib = which( as.numeric( rownames( B) ) %in% spring )
      B[ib,3] = NA 

      L = L[ which(rownames(L) %in% rownames(B) ),  ]
      B.sd = B.sd[ which(rownames(B.sd) %in% rownames(B) ),  ]
      R = R[ which(rownames(R) %in% rownames(B) ),  ]
      R.sd = R.sd[ which(rownames(R.sd) %in% rownames(B) ),  ]
      Bg = Bg[ which(rownames(Bg) %in% rownames(B) ),  ]
      Bg.sd = Bg.sd[ which(rownames(Bg.sd) %in% rownames(B) ),  ]
      Rg = Rg[ which(rownames(Rg) %in% rownames(B) ),  ]
      Rg.sd = Rg.sd[ which(rownames(Rg.sd) %in% rownames(B) ),  ]
  
      out = list( L=L, B=B, R=R, B.sd=B.sd, R.sd=R.sd, Bg=Bg, Bg.sd=Bg.sd, Rg=Rg, Rg.sd=Rg.sd ) 

      save( out, file=fn, compress=T )
      
      return( fn )
    }


    if (DS %in% c("L", "L.redo" ) ) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "L.bugs.rdata" )
      L=NULL
      if (DS=="L") {
        if (file.exists(fn)) load(fn)
        return(L)
      } 
      L = landings.aggregate( format="bugs" )
      L = L[ , c("cfanorth", "cfasouth", "cfa4x") ]
      L = as.data.frame(L)
      save( L, file=fn, compress=T )
      return (L)
    }


    if (DS %in% c("B", "B.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "B.bugs.rdata" ) 
      B=NULL
      if (DS=="B") {
        if (file.exists(fn)) load(fn)
        return(B)
      } 
       
      # biomass data: post-fishery biomass are determined by survey B)
        p$vars.to.model ="R0.mass"
        p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
        
        K = interpolation.db( DS="interpolation.simulation", p=p )
        areas=c("cfanorth", "cfasouth", "cfa4x")
        td = K[ which( K$region %in% areas) ,]

        B = tapply( td$total, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        B = B / 1000 # kt
        B = B[ , areas]
        B = as.data.frame(B)
        save( B, file=fn, compress=T )
        return (B)
    
    }
   
    
    if (DS %in% c("B.sd", "B.sd.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "B.sd.bugs.rdata" ) 
      B=NULL
      if (DS=="B.sd") {
        if (file.exists(fn)) load(fn)
        return(B)
      } 
       
      # biomass data: post-fishery biomass are determined by survey B)
        p$vars.to.model ="R0.mass"
        p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
        
        K = interpolation.db( DS="interpolation.simulation", p=p )
        areas=c("cfanorth", "cfasouth", "cfa4x")
        td = K[ which( K$region %in% areas) ,]

        B = tapply( td$total.sd.ln, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        B = B[ , areas]
        B = as.data.frame(B)
        save( B, file=fn, compress=T )
        return (B)
    }



    if (DS %in% c("R_geomean", "R_geomean.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "R_geomean.bugs.rdata" ) 
      Bx = NULL
      if (DS=="R_geomean") {
        if (file.exists(fn)) load(fn)
        return(Bx)
      } 

      v = "R1.no" 
      fm = formula( paste(v, "~yr+cfa"))
      
      set = snowcrab.db( DS ="set.complete", p=p, yrs=1996:p$current.assessment.year ) 
      B = set[, c( v, "yr", "cfa" )]
      B[,v] = variable.recode( B[,v], v, direction="forward", db="snowcrab" )

      Bg = aggregate( fm , data=B, FUN=mean, na.rm=T )
      Bg[,v] = variable.recode( Bg[,v], v, direction="backward", db="snowcrab" )

      Bx = xtabs( fm, data=Bg )

      Bx[ which( rownames(Bx) %in% c(1996:2001) ),"cfa4x"] = NA
      Bx[ which( rownames(Bx) %in% c(1996) ),"cfanorth"] = NA
      Bx = Bx[,c("cfanorth", "cfasouth", "cfa4x" ) ] 
      
      Bx[ which(Bx < 0.0001 )] = 0
 
      save( Bx, file=fn, compress=T )
      return (Bx)
    }
 


    if (DS %in% c("B_geomean", "B_geomean.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "B_geomean.bugs.rdata" ) 
      Bx = NULL
      if (DS=="B_geomean") {
        if (file.exists(fn)) load(fn)
        return(Bx)
      } 

      v = "R0.mass" 
      fm = formula( paste(v, "~yr+cfa"))
      
      set = snowcrab.db( DS ="set.complete", p=p, yrs=1996:p$current.assessment.year ) 
      B = set[, c( v, "yr", "cfa" )]
      B[,v] = variable.recode( B[,v], v, direction="forward", db="snowcrab" )

      Bg = aggregate( fm , data=B, FUN=mean, na.rm=T )
      Bg[,v] = variable.recode( Bg[,v], v, direction="backward", db="snowcrab" )

      Bx = xtabs( fm, data=Bg )

      Bx[ which( rownames(Bx) %in% c(1996:2001) ),"cfa4x"] = NA
      Bx[ which( rownames(Bx) %in% c(1996) ),"cfanorth"] = NA
      Bx = Bx[,c("cfanorth", "cfasouth", "cfa4x" ) ] 

      save( Bx, file=fn, compress=T )
      return (Bx)
    }
 

    if (DS %in% c("B_geomean.sd", "B_geomean.sd.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "B_geomean.sd.bugs.rdata" ) 
      Bx = NULL
      if (DS=="B_geomean.sd") {
        if (file.exists(fn)) load(fn)
        return(Bx)
      } 

      v = "R0.mass" 
      fm = formula( paste(v, "~yr+cfa"))
      
      set = snowcrab.db( DS ="set.complete", p=p, yrs=1996:p$current.assessment.year ) 
      B = set[, c( v, "yr", "cfa" )]
      B[,v] = variable.recode( B[,v], v, direction="forward", db="snowcrab" )

      Bg = aggregate( fm , data=B, FUN=sd, na.rm=T )
      # Bg[,v] = variable.recode( Bg[,v], v, direction="backward", db="snowcrab" )

      Bx = xtabs( fm, data=Bg )

      Bx[ which( rownames(Bx) %in% c(1996:2001) ),"cfa4x"] = NA
      Bx[ which( rownames(Bx) %in% c(1996) ),"cfanorth"] = NA
      Bx = Bx[, c("cfanorth", "cfasouth", "cfa4x" ) ] 

      save( Bx, file=fn, compress=T )
      return (Bx)
    }
 

    if (DS %in% c("R_geomean.sd", "R_geomean.sd.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "R_geomean.sd.bugs.rdata" ) 
      Bx = NULL
      if (DS=="R_geomean.sd") {
        if (file.exists(fn)) load(fn)
        return(Bx)
      } 

      v = "R1.no" 
      fm = formula( paste(v, "~yr+cfa"))
      
      set = snowcrab.db( DS ="set.complete", p=p, yrs=1996:p$current.assessment.year ) 
      B = set[, c( v, "yr", "cfa" )]
      B[,v] = variable.recode( B[,v], v, direction="forward", db="snowcrab" )

      Bg = aggregate( fm , data=B, FUN=sd, na.rm=T )
      # Bg[,v] = variable.recode( Bg[,v], v, direction="backward", db="snowcrab" )

      Bx = xtabs( fm, data=Bg )

      Bx[ which( rownames(Bx) %in% c(1996:2001) ),"cfa4x"] = NA
      Bx[ which( rownames(Bx) %in% c(1996) ),"cfanorth"] = NA
      Bx = Bx[,c("cfanorth", "cfasouth", "cfa4x" ) ] 

      save( Bx, file=fn, compress=T )
      return (Bx)
    }
 
    if (DS %in% c("R", "R.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "R.bugs.rdata" ) 
      R=NULL
      if (DS=="R") {
        if (file.exists(fn)) load(fn)
        return(R)
      } 
       
        p$vars.to.model = "R1.no"
        p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
        
        K = interpolation.db( DS="interpolation.simulation", p=p )
        areas=c("cfanorth", "cfasouth", "cfa4x")
        td = K[ which(K$vars==p$vars.to.model & K$region %in% areas) ,]
        R = tapply( td$total, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        R = R /1000 # kt
        R = R[ , areas]

      # assessments conducted in the spring 2001 and earlier ... they are pre-fishery biomasses .. 
      # fishery is modelled using a biomass (fishable biomass) estimate == pre-fishery biomass == landings + post-fishery biomass
        
        R = R[, areas]
        R[ which( rownames(R) %in% c(1998:2001) ),"cfa4x"] = mean( R[ which( rownames(R) %in% c(2002:2005) ),"cfa4x"] )
        save( R, file=fn, compress=T )
        return (R)
    }



    if (DS %in% c("R.sd", "R.sd.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "R.sd.bugs.rdata" ) 
      R=NULL
      if (DS=="R.sd") {
        if (file.exists(fn)) load(fn)
        return(R)
      } 
       
        p$vars.to.model = "R1.no"
        p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
        
        K = interpolation.db( DS="interpolation.simulation", p=p )
        areas=c("cfanorth", "cfasouth", "cfa4x")
        td = K[ which(K$vars==p$vars.to.model & K$region %in% areas) ,]
        R = tapply( td$total.sd.ln, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        R = R[ , areas]

      # assessments conducted in the spring 2001 and earlier ... they are pre-fishery biomasses .. 
      # fishery is modelled using a biomass (fishable biomass) estimate == pre-fishery biomass == landings + post-fishery biomass
        
        R[ which( rownames(R) %in% c(1998:2001) ),"cfa4x"] = mean( R[ which( rownames(R) %in% c(2002:2005) ),"cfa4x"] )
        R[ which( R < 0.001) ] = mean( R, na.rm=T )
        
        save( R, file=fn, compress=T )
        return (R)
    }

  }



