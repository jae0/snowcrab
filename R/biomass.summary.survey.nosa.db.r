
  biomass.summary.survey.nosa.db = function( DS="complete", p=NULL ) {
    #browser()
    
    if (DS %in% c("complete", "complete.redo") ) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "bugssurvey.rdata" ) 
      if ( DS == "complete") {
        out = NULL
        if (file.exists(fn)) load(fn)
        return( out )
      }
      L = biomass.summary.survey.nosa.db( DS="L.redo", p=p  )  # must go first as part of biomass estimates
      B = biomass.summary.survey.nosa.db( DS="B.redo", p=p )  # rename to avoid confusion below as B is also used
      B.sd = biomass.summary.survey.nosa.db( DS="B.sd.redo", p=p )  # rename to avoid confusion below as B is also used
      R = biomass.summary.survey.nosa.db( DS="R.redo", p=p  )  # rename to avoid confusion below as B is also used
      R.sd = biomass.summary.survey.nosa.db( DS="R.sd.redo", p=p  )  # rename to avoid confusion below as B is also used
      
      # geometric means -- sd are on log scale
      Bg = biomass.summary.survey.nosa.db( DS="B_geomean.redo", p=p  )
      Bg.sd = biomass.summary.survey.nosa.db( DS="B_geomean.sd.redo", p=p  )
      Rg = biomass.summary.survey.nosa.db( DS="R_geomean.redo", p=p  )
      Rg.sd = biomass.summary.survey.nosa.db( DS="R_geomean.sd.redo", p=p  )
 

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
      fn = file.path( project.datadirectory("snowcrab"), "R", "L.bugs.survey.rdata" )
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
      fn = file.path( project.datadirectory("snowcrab"), "R", "B.bugs.survey.rdata" ) 
      B=NULL
      if (DS=="B") {
        if (file.exists(fn)) load(fn)
        return(B)
      } 
       
      # biomass data: post-fishery biomass are determined by survey B)
        p$vars.to.model ="R0.mass"
        p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
        
        #fdir <-file.path( project.datadirectory("snowcrab"), "R", "gam","habitat" )
        #fs <- dir(fdir)
		    #fs <- fs[setdiff(grep('K.R0.mass',fs) , grep('environmentals.only',fs))]
        #lo <- c()
        #for(i in 1:length(fs)) {
        #	load(file.path(fdir,fs[i]))
        #lo <- rbind(lo,K)
        #rm(K)
        #}
        
        areas=c("cfanorth", "cfasouth", "cfa4x")
        #td = lo[ which( lo$region %in% areas) , c('yr','region','sa.region')]
        
	  	load(file.path(project.datadirectory('snowcrab'),"R","ts.rdata"))
		  ts <- ts[which(ts$variable==p$vars.to.model & ts$region %in% areas),c('year','region','mean')]
		  ts$yr <- ts$year
		#td <- merge(ts,td, by=c('yr','region'),all.x=T)
		#if(any(is.na(td$sa.region))) {
		#j <- td[which(is.na(td$sa.region)),]
		#id <- unique(j[,c('region','yr')])
		#	for(i in 1:nrow(id)) {
		#	yy <- id[i,'yr']
		#	y <- (yy-5):(yy-1)
		#	td[which(td$region==id[i,'region'] & td$yr == id[i,'yr']),'sa.region']	 <- mean(td[which(td$region==id[i,'region'] & td$yr %in% y),'sa.region'],na.rm=T)
		#	}
		#}
	#	td$total <- td$mean*td$sa.region
   #     B = tapply( td$total, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
    #    B = B / 1000 # kt
     #   B = B[ , areas]
        
        B = tapply( ts$mean, INDEX=ts[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        #B = B / 1000 # kt
        B = B[ , areas]
        B = as.data.frame(B)
        print('Biomass Index Survey')
        print(B)
        save( B, file=fn, compress=T )
        return (B)
       
    }
   
    
    if (DS %in% c("B.sd", "B.sd.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "B.sd.bugs.survey.rdata" ) 
      B=NULL
      if (DS=="B.sd") {
        if (file.exists(fn)) load(fn)
        return(B)
      } 
       
        p$vars.to.model ="R0.mass"
        p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
        
        #fdir <-file.path( project.datadirectory("snowcrab"), "R", "gam","habitat" )
        #fs <- dir(fdir)
		    #fs <- fs[setdiff(grep('K.R0.mass',fs) , grep('environmentals.only',fs))]
        #lo <- c()
        #for(i in 1:length(fs)) {
        #	load(file.path(fdir,fs[i]))
        #lo <- rbind(lo,K)
        #rm(K)
        #}
        
        areas=c("cfanorth", "cfasouth", "cfa4x")
        #td = lo[ which( lo$region %in% areas) , c('yr','region','sa.region')]
        
		load(file.path(project.datadirectory('snowcrab'),"R","ts.rdata"))
		ts <- ts[which(ts$variable==p$vars.to.model & ts$region %in% areas),c('year','region','se')]
		ts$yr <- ts$year
		#td <- merge(ts,td, by=c('yr','region'),all.x=T)
		#if(any(is.na(td$sa.region))) {
		#j <- td[which(is.na(td$sa.region)),]
		#id <- unique(j[,c('region','yr')])
		#	for(i in 1:nrow(id)) {
		#	yy <- id[i,'yr']
		#	y <- (yy-5):(yy-1)
		#	td[which(td$region==id[i,'region'] & td$yr == id[i,'yr']),'sa.region']	 <- mean(td[which(td$region==id[i,'region'] & td$yr %in% y),'sa.region'],na.rm=T)
		#	}
		#}
		#td$total <- td$se*td$sa.region
        B = tapply( ts$se, INDEX=ts[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        B = log(B+1) # kt
        B = B[ , areas]
        B = as.data.frame(B)
        save( B, file=fn, compress=T )
        return (B)
    }



    if (DS %in% c("R_geomean", "R_geomean.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "R_geomean.bugs.survey.rdata" ) 
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
      fn = file.path( project.datadirectory("snowcrab"), "R", "B_geomean.bugs.survey.rdata" ) 
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
      fn = file.path( project.datadirectory("snowcrab"), "R", "B_geomean.sd.bugs.survey.rdata" ) 
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
      fn = file.path( project.datadirectory("snowcrab"), "R", "R_geomean.sd.bugs.survey.rdata" ) 
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
      fn = file.path( project.datadirectory("snowcrab"), "R", "R.bugs.survey.rdata" ) 
      R=NULL
      if (DS=="R") {
        if (file.exists(fn)) load(fn)
        return(R)
      } 
       
        p$vars.to.model = "R1.no"
        p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
                p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
        
        fdir <-file.path( project.datadirectory("snowcrab"), "R", "gam","habitat" )
        fs <- dir(fdir)
		fs <- fs[setdiff(grep('K.R0.mass',fs) , grep('environmentals.only',fs))]
        lo <- c()
        for(i in 1:length(fs)) {
        	load(file.path(fdir,fs[i]))
        lo <- rbind(lo,K)
        rm(K)
        }
        
        areas=c("cfanorth", "cfasouth", "cfa4x")
        td = lo[ which( lo$region %in% areas) , c('yr','region','sa.region')]
        
		load(file.path(project.datadirectory('snowcrab'),"R","ts.rdata"))
		ts <- ts[which(ts$variable==p$vars.to.model & ts$region %in% areas),c('year','region','mean')]
    ts[ts$mean<0.1,'mean'] <- min(ts$mean[ts$mean>0.1])
		ts$yr <- ts$year
		td <- merge(ts,td, by=c('yr','region'),all.x=T)
		if(any(is.na(td$sa.region))) {
		j <- td[which(is.na(td$sa.region)),]
		id <- unique(j[,c('region','yr')])
			for(i in 1:nrow(id)) {
			yy <- id[i,'yr']
			y <- (yy-5):(yy-1)
			td[which(td$region==id[i,'region'] & td$yr == id[i,'yr']),'sa.region']	 <- mean(td[which(td$region==id[i,'region'] & td$yr %in% y),'sa.region'],na.rm=T)
			}
		}
		td$total <- td$mean*td$sa.region
        R = tapply( td$total, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        R = R / 1000 # kt
        R = R[ , areas]
        
       
        
        R[ which( rownames(R) %in% c(1998:2001) ),"cfa4x"] = mean( R[ which( rownames(R) %in% c(2002:2005) ),"cfa4x"] )
        save( R, file=fn, compress=T )
      # assessments conducted in the spring 2001 and earlier ... they are pre-fishery biomasses .. 
      # fishery is modelled using a biomass (fishable biomass) estimate == pre-fishery biomass == landings + post-fishery biomass
      return(R)
      }



    if (DS %in% c("R.sd", "R.sd.redo" )) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "R.sd.bugs.survey.rdata" ) 
      R=NULL
      if (DS=="R.sd") {
        if (file.exists(fn)) load(fn)
        return(R)
      } 
       
               p$vars.to.model = "R1.no"
        p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
                p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
        
        fdir <-file.path( project.datadirectory("snowcrab"), "R", "gam","habitat" )
        fs <- dir(fdir)
		fs <- fs[setdiff(grep('K.R0.mass',fs) , grep('environmentals.only',fs))]
        lo <- c()
        for(i in 1:length(fs)) {
        	load(file.path(fdir,fs[i]))
        lo <- rbind(lo,K)
        rm(K)
        }
        
        areas=c("cfanorth", "cfasouth", "cfa4x")
        td = lo[ which( lo$region %in% areas) , c('yr','region','sa.region')]
        
		load(file.path(project.datadirectory('snowcrab'),"R","ts.rdata"))
		ts <- ts[which(ts$variable==p$vars.to.model & ts$region %in% areas),c('year','region','se')]
		ts$yr <- ts$year
		td <- merge(ts,td, by=c('yr','region'),all.x=T)
		if(any(is.na(td$sa.region))) {
		j <- td[which(is.na(td$sa.region)),]
		id <- unique(j[,c('region','yr')])
			for(i in 1:nrow(id)) {
			yy <- id[i,'yr']
			y <- (yy-5):(yy-1)
			td[which(td$region==id[i,'region'] & td$yr == id[i,'yr']),'sa.region']	 <- mean(td[which(td$region==id[i,'region'] & td$yr %in% y),'sa.region'],na.rm=T)
			}
		}
		td$total <- td$se*td$sa.region
        R = tapply( td$total, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        
        R = R[ , areas]
        #R = as.data.frame(R)
        R[ which( rownames(R) %in% c(1996:2001) ),"cfa4x"] = mean( R[ which( rownames(R) %in% c(2002:2005) ),"cfa4x"] )
        R[ which( rownames(R) %in% c(1996) ),"cfanorth"] = mean( R[ which( rownames(R) %in% c(1997:2002 )),"cfanorth"] )
       
         R[ which( R < 0.001) ] = sum(R)/(dim(R)[1]*dim(R)[2])
      
        save( R, file=fn, compress=T )
      # assessments conducted in the spring 2001 and earlier ... they are pre-fishery biomasses .. 
      # fishery is modelled using a biomass (fishable biomass) estimate == pre-fishery biomass == landings + post-fishery biomass
      return(R)
        
    }

  }




  