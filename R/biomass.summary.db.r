
  biomass.summary.db = function( DS="complete", p=NULL ) {
    #browser()
    
    sum_outdir = file.path( project.datadirectory("bio.snowcrab"), "output", "timeseries" )

    if (DS=="surplusproduction" ){
      res =  biomass.summary.db(p=p, DS="complete")
      sb = list(
        bmin = 0.001, # scaled to 1 but allow overshooting
        bmax = 1.1, # scaled to 1 but allow overshooting
        qmin = 0.001,
        qmax = 2,
        qmu = rep(1,3),
        qsd = rep(0.3,3),
        bomup=rep(-2.5579,3),
        bosdp=rep(0.47726,3),
        bpmup=rep(-2.5579,3),
        bpsdp=rep(0.47726,3),
        bomin=rep(0,3),
        bomax=rep(5,3),
        bpmin=rep(0,3),
        bpmax=rep(5,3),
        recmax= c( 10^3, 10^4, 10^2 ),
        Kmu = c( 1.831139,4.170013,0.784308), #for ln
        Ksd = c(0.06,0.06,0.04), #for ln
        rmu = rep(0.96,3),
        rsd = rep(0.01041271,3),
        b0min = c(0.5, 0.5, 0.2),  # prior: mean value possible in  N,S,4X
        b0max = c(0.8, 0.8, 0.6),  # prior: mean value possible in  N,S,4X
        cvnormalmin = 0.05, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
        cvnormalmax = 0.4, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
        cvlognormalmin = 0.05, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
        cvlognormalmax = 0.4, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
      # for lognormal: cv = sqrt(exp(sigma^2) - 1); or sigma = sqrt(log(cv^2+ 1) ) ==> sigma = sqrt( log(0.25^2 + 1)) = 0.246 ~ cv -- i.e. cv ~ sd
        IOA = as.matrix(res$B), # observed index of abundance
        IOAcv = as.matrix( log(res$B.sd )), # observed index of log abundance SD estimates ~ CV
        #IREC = as.matrix(res$R), # observed index of abundance
        #IRECcv = as.matrix(res$Rsd ), # observed index of log abundance SD estimates ~CV
        CAT = as.matrix(res$L) , # catches  , assume 20% handling mortality and illegal landings
        CATmin = apply( res$L, 2, min, na.rm=T),
        CATmax = apply( res$L, 2, max, na.rm=T),
        er = 0.2,  # target exploitation rate
        U = ncol( res$B),  # number of regions
        N = nrow( res$B) , # no years with data
        M = 3, # no years for projections
        ty = 7,  # index of the transition year (2004) between spring and fall surveys
        cfa4x = 3, # index of cfa4x
        eps = 1e-4  # small non-zero number
      )

    # set poor data to NA's and maximize associated CV's

      cfa.north =  1 # column index
      cfa.north.bad.data = which( as.numeric(rownames(sb$IOA)) <= 1997 )
      sb$IOA[ cfa.north.bad.data, cfa.north ] = NA
      sb$IOAcv[ cfa.north.bad.data, cfa.north ] = mean( sb$IOAcv[,cfa.north ] ,na.rm=T )
      #sb$IRECcv[ cfa.north.bad.data, cfa.north ] = mean( sb$IRECcv[,cfa.north ] ,na.rm=T )
      #sb$IREC[ cfa.north.bad.data, cfa.north ] = mean( sb$IREC[,cfa.north ] ,na.rm=T )

      cfa.south =  2 # column index
      cfa.south.bad.data = which( as.numeric(rownames(sb$IOA)) <= 1998 )
      sb$IOA[ cfa.south.bad.data, cfa.south ] = NA
      sb$IOAcv[ cfa.south.bad.data, cfa.south ] = mean( sb$IOAcv[,cfa.south ] ,na.rm=T )

      cfa.nodata =   which( as.numeric(rownames(sb$IOA)) <= 2003 )
      sb$IOA[ cfa.nodata , sb$cfa4x ] = NA
      sb$IOAcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IOAcv[,sb$cfa4x] ,na.rm=T )

      # recruitment index has no process model right now ... use min value
      #sb$IREC[ cfa.nodata , sb$cfa4x ] = min( sb$IREC[,sb$cfa4x] ,na.rm=T )
      #sb$IRECcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IRECcv[,sb$cfa4x] ,na.rm=T )

      #sb$IREC = log( sb$IREC )
      return(sb)
    }


    if (DS %in% c("complete", "complete.redo") ) {
      fn = file.path( sum_outdir, "complete_ts.rdata" )
      if ( DS == "complete") {
        out = NULL
        if (file.exists(fn)) load(fn)
        return( out )
      }

      L = biomass.summary.db( DS="L.redo", p=p  )  # must go first as part of biomass estimates
      B = biomass.summary.db( DS="B.redo", p=p )  # rename to avoid confusion below as B is also used
      Bsd = biomass.summary.db( DS="Bsd.redo", p=p )  # rename to avoid confusion below as B is also used
     #  R = biomass.summary.db( DS="R.redo", p=p  )  # rename to avoid confusion below as B is also used
     #  Rsd = biomass.summary.db( DS="Rsd.redo", p=p  )  # rename to avoid confusion below as B is also used

      # geometric means -- sd are on log scale
      Bg = biomass.summary.db( DS="B_geomean.redo", p=p  )
      Bgsd = biomass.summary.db( DS="B_geomeansd.redo", p=p  )
     # Rg = biomass.summary.db( DS="R_geomean.redo", p=p  )
     # Rgsd = biomass.summary.db( DS="R_geomeansd.redo", p=p  )


      # cfa4x have had no estimates prior to 2004
      spring = 1998:2003
      ib = which( as.numeric( rownames( B) ) %in% spring )
      B[ib,3] = NA

      L = L[ which(rownames(L) %in% rownames(B) ),  ]
      Bsd = Bsd[ which(rownames(Bsd) %in% rownames(B) ),  ]
    #  R = R[ which(rownames(R) %in% rownames(B) ),  ]
    #  Rsd = Rsd[ which(rownames(Rsd) %in% rownames(B) ),  ]
      Bg = Bg[ which(rownames(Bg) %in% rownames(B) ),  ]
      Bgsd = Bgsd[ which(rownames(Bgsd) %in% rownames(B) ),  ]
    #  Rg = Rg[ which(rownames(Rg) %in% rownames(B) ),  ]
    #  Rgsd = Rgsd[ which(rownames(Rgsd) %in% rownames(B) ),  ]

     out = list( L=L, B=B, Bsd=Bsd,Bg=Bg, Bgsd=Bgsd )

    #  out = list( L=L, B=B, R=R, Bsd=Bsd, Rsd=Rsd, Bg=Bg, Bgsd=Bgsd, Rg=Rg, Rgsd=Rgsd )

      save( out, file=fn, compress=T )

      return( fn )
    }


    if (DS %in% c("L", "L.redo" ) ) {
      fn = file.path( sum_outdir, "L_ts.rdata" )
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
      fn = file.path( sum_outdir, "B_ts.rdata" )
      B=NULL
      if (DS=="B") {
        if (file.exists(fn)) load(fn)
        return(B)
      }

      # biomass data: post-fishery biomass are determined by survey B)
        
        K = interpolation.db( DS="interpolation.simulation", p=p, varnames="R0.mass" )
        areas=c("cfanorth", "cfasouth", "cfa4x")
        td = K[ which( K$region %in% areas) ,]

        B = tapply( td$total, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        B = B / 1000 # kt
        B = B[ , areas]
        B = as.data.frame(B)
        save( B, file=fn, compress=T )
        return (B)

    }


    if (DS %in% c("Bsd", "Bsd.redo" )) {
      fn = file.path( sum_outdir, "Bsd_ts.rdata" )
      B=NULL
      if (DS=="Bsd") {
        if (file.exists(fn)) load(fn)
        return(B)
      }
      # biomass data: post-fishery biomass are determined by survey B)
        K = interpolation.db( DS="interpolation.simulation", p=p )
        areas=c("cfanorth", "cfasouth", "cfa4x")
        td = K[ which( K$region %in% areas) ,]
        td$totalsd.ln = abs(td$total.ub - td$total.lb)/(2*1.96) # assuming approx normal
        B = tapply( td$totalsd.ln, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
        B = B[ , areas]
        B = as.data.frame(B)
        save( B, file=fn, compress=T )
        return (B)
    }



    if (DS %in% c("R_geomean", "R_geomean.redo" )) {
      fn = file.path( sum_outdir, "R_geomean_ts.rdata" )
      Bx = NULL
      if (DS=="R_geomean") {
        if (file.exists(fn)) load(fn)
        return(Bx)
      }

      lookup.table = snowcrab.db( p=p, DS="data.transforms" )

      v = "R1.no"
      fm = formula( paste(v, "~yr+cfa"))

      set = snowcrab.db( DS ="set.biologicals", p=p, yrs=1996:p$year.assessment )
      B = set[, c( v, "yr", "cfa" )]
      B[,v] = bio.snowcrab::variable.recode( B[,v], v, direction="forward", lookup.table=lookup.table )

      Bg = aggregate( fm , data=B, FUN=mean, na.rm=T )
      Bg[,v] = bio.snowcrab::variable.recode( Bg[,v], v, direction="backward", lookup.table=lookup.table )

      Bx = xtabs( fm, data=Bg )

      Bx[ which( rownames(Bx) %in% c(1996:2001) ),"cfa4x"] = NA
      Bx[ which( rownames(Bx) %in% c(1996) ),"cfanorth"] = NA
      Bx = Bx[,c("cfanorth", "cfasouth", "cfa4x" ) ]

      Bx[ which(Bx < 0.0001 )] = 0

      save( Bx, file=fn, compress=T )
      return (Bx)
    }



    if (DS %in% c("B_geomean", "B_geomean.redo" )) {
      fn = file.path( sum_outdir, "B_geomean_ts.rdata" )
      Bx = NULL
      if (DS=="B_geomean") {
        if (file.exists(fn)) load(fn)
        return(Bx)
      }
      
      lookup.table = snowcrab.db( p=p, DS="data.transforms" )

      v = "R0.mass"
      fm = formula( paste(v, "~yr+cfa"))

      set = snowcrab.db( DS ="set.biologicals", p=p, yrs=1996:p$year.assessment )
      B = set[, c( v, "yr", "cfa" )]
      B[,v] = bio.snowcrab::variable.recode( B[,v], v, direction="forward", lookup.table=lookup.table )

      Bg = aggregate( fm , data=B, FUN=mean, na.rm=T )
      Bg[,v] = bio.snowcrab::variable.recode( Bg[,v], v, direction="backward", lookup.table=lookup.table )

      Bx = xtabs( fm, data=Bg )

      Bx[ which( rownames(Bx) %in% c(1996:2001) ),"cfa4x"] = NA
      Bx[ which( rownames(Bx) %in% c(1996) ),"cfanorth"] = NA
      Bx = Bx[,c("cfanorth", "cfasouth", "cfa4x" ) ]

      save( Bx, file=fn, compress=T )
      return (Bx)
    }


    if (DS %in% c("B_geomeansd", "B_geomeansd.redo" )) {
      fn = file.path( sum_outdir, "B_geomeansd_ts.rdata" )
      Bx = NULL
      if (DS=="B_geomeansd") {
        if (file.exists(fn)) load(fn)
        return(Bx)
      }
      
      lookup.table = snowcrab.db( p=p, DS="data.transforms" )

      v = "R0.mass"
      fm = formula( paste(v, "~yr+cfa"))

      set = snowcrab.db( DS ="set.biologicals", p=p, yrs=1996:p$year.assessment )
      B = set[, c( v, "yr", "cfa" )]
      B[,v] = bio.snowcrab::variable.recode( B[,v], v, direction="forward", lookup.table=lookup.table )

      Bg = aggregate( fm , data=B, FUN=sd, na.rm=T )
      # Bg[,v] = bio.snowcrab::variable.recode( Bg[,v], v, direction="backward", lookup.table=lookup.table )

      Bx = xtabs( fm, data=Bg )

      Bx[ which( rownames(Bx) %in% c(1996:2001) ),"cfa4x"] = NA
      Bx[ which( rownames(Bx) %in% c(1996) ),"cfanorth"] = NA
      Bx = Bx[, c("cfanorth", "cfasouth", "cfa4x" ) ]

      save( Bx, file=fn, compress=T )
      return (Bx)
    }


    if (DS %in% c("R_geomeansd", "R_geomeansd.redo" )) {
      fn = file.path( sum_outdir, "R_geomeansd_ts.rdata" )
      Bx = NULL
      if (DS=="R_geomeansd") {
        if (file.exists(fn)) load(fn)
        return(Bx)
      }
     
      lookup.table = snowcrab.db( p=p, DS="data.transforms" )

      v = "R1.no"
      fm = formula( paste(v, "~yr+cfa"))

      set = snowcrab.db( DS ="set.biologicals", p=p, yrs=1996:p$year.assessment )
      B = set[, c( v, "yr", "cfa" )]
      B[,v] = bio.snowcrab::variable.recode( B[,v], v, direction="forward", lookup.table=lookup.table )

      Bg = aggregate( fm , data=B, FUN=sd, na.rm=T )
      # Bg[,v] = bio.snowcrab::variable.recode( Bg[,v], v, direction="backward", lookup.table=lookup.table )

      Bx = xtabs( fm, data=Bg )

      Bx[ which( rownames(Bx) %in% c(1996:2001) ),"cfa4x"] = NA
      Bx[ which( rownames(Bx) %in% c(1996) ),"cfanorth"] = NA
      Bx = Bx[,c("cfanorth", "cfasouth", "cfa4x" ) ]

      save( Bx, file=fn, compress=T )
      return (Bx)
    }

    if (DS %in% c("R", "R.redo" )) {
      fn = file.path( sum_outdir, "R_ts.rdata" )
      R=NULL
      if (DS=="R") {
        if (file.exists(fn)) load(fn)
        return(R)
      }

      K = interpolation.db( DS="interpolation.simulation", p=p, varnames="R1.no" )
      areas=c("cfanorth", "cfasouth", "cfa4x")
      td = K[ which(K$vars=="R1.no" & K$region %in% areas) ,]
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



    if (DS %in% c("Rsd", "Rsd.redo" )) {
      fn = file.path( sum_outdir, "Rsd_ts.rdata" )
      R=NULL
      if (DS=="Rsd") {
        if (file.exists(fn)) load(fn)
        return(R)
      }
      message ("not available right now")

      stop()
      
      K = interpolation.db( DS="interpolation.simulation", p=p, varnames="R1.no" )
      areas=c("cfanorth", "cfasouth", "cfa4x")
      td = K[ which(K$vars=="R1.no" & K$region %in% areas) ,]
      td$totalsd.ln = abs(td$total.ub - td$total.lb)/(2*1.96) # assuming approx normal
      R = tapply( td$totalsd.ln, INDEX=td[,c("yr", "region")], FUN=sum, na.rm=T )  # summation is really returning identity as there is only 1 element
      R = R[ , areas]

    # assessments conducted in the spring 2001 and earlier ... they are pre-fishery biomasses ..
    # fishery is modelled using a biomass (fishable biomass) estimate == pre-fishery biomass == landings + post-fishery biomass

      R[ which( rownames(R) %in% c(1998:2001) ),"cfa4x"] = mean( R[ which( rownames(R) %in% c(2002:2005) ),"cfa4x"] )
      R[ which( R < 0.001) ] = mean( R, na.rm=T )

      save( R, file=fn, compress=T )
      return (R)
    }

  }



