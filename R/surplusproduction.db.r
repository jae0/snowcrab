
surplusproduction.db = function( DS, sourcedata="default", debug.region="cfanorth" ) {
    # sb = surplusproduction.db( DS="jags.2015", sourcedata="nosa" ) 
    if (is.data.frame(sourcedata)) {
      res= sourcedata 
    } else {
      if (sourcedata=="default") res =  biomass.summary.db(p=p)
      if (sourcedata=="survey") res =  biomass.summary.survey.db(p=p)
      if (sourcedata=="nosa" )   res = biomass.summary.survey.nosa.db(p=p)
    }

  #\\ Create data required for surplus production modelling
  #\\ NOTE: for lognormal: cv = sqrt(exp(sigma^2) - 1); 
  #\\ or sigma = sqrt(log(cv^2+ 1) ) ==> sigma = sqrt( log(0.25^2 + 1)) = 0.246 ~ cv -- i.e. cv ~ sd
  if (DS %in% c("jags.2013")) {

    sb = list(
      b.min = 0.001, # scaled to 1 but allow overshooting
      b.max = 1.1, # scaled to 1 but allow overshooting
      q.min = 0.1,  # max value of q , anything larger is not believable
      q.max = 2.0,  # max value of q , anything larger is not believable
      r.min = 0.1,
      r.max = 2.0,
      rec.max= c( 10^3, 10^4, 10^2 ),
      K.min = c(1, 10, 0.1 ),  # max carrying capacity estimate:
      K.max = c(10, 100, 5 ),  # max carrying capacity estimate:
      b0.min = c(0.5, 0.5, 0.2),  # prior: mean value possible in  N,S,4X
      b0.max = c(0.8, 0.8, 0.6),  # prior: mean value possible in  N,S,4X
      cv.normal.min = 0.05, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
      cv.normal.max = 0.4, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
      cv.lognormal.min = 0.05, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
      cv.lognormal.max = 0.4, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
      IOA = as.matrix(res$B), # observed index of abundance
      IOAcv = as.matrix(res$B.sd ), # observed index of log abundance SD estimates ~ CV
      IREC = as.matrix(res$R), # observed index of abundance
      IRECcv = as.matrix(res$R.sd ), # observed index of log abundance SD estimates ~CV
      CAT = as.matrix(res$L) , # catches  , assume 20% handling mortality and illegal landings
      CAT.min = apply( res$L, 2, min, na.rm=T),
      CAT.max = apply( res$L, 2, max, na.rm=T),
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

    cfa.south =  2 # column index
    cfa.south.bad.data = which( as.numeric(rownames(sb$IOA)) <= 1998 )
    sb$IOA[ cfa.south.bad.data, cfa.south ] = NA
    sb$IOAcv[ cfa.south.bad.data, cfa.south ] = mean( sb$IOAcv[,cfa.south ] ,na.rm=T )

    cfa.nodata =   which( as.numeric(rownames(sb$IOA)) <= 2003 )
    sb$IOA[ cfa.nodata , sb$cfa4x ] = NA
    sb$IOAcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IOAcv[,sb$cfa4x] ,na.rm=T )

    # recruitment index has no process model right now ... use min value
    sb$IREC[ cfa.nodata , sb$cfa4x ] = min( sb$IREC[,sb$cfa4x] ,na.rm=T )
    sb$IRECcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IRECcv[,sb$cfa4x] ,na.rm=T )

    sb$IREC = log( sb$IREC )

    sb$tomonitor = c( "r", "K", "q", "qs", "r.mu", "r.sd", "b","bp.sd", "bo.sd", "b0", "b0.sd", "F", "TAC",  "C", "P", "B", "rem", "rem.sd", "rem.mu", "REM", "MSY", "BMSY", "FMSY", "Fcrash", "Bdrop", "BX2MSY" ) 
    sb$jagsmodelname = "biomassdynamic_full_2013.bugs"

    return(sb)
  }



  if ( DS %in% c("jags.2014")  ) {

    sb = list(
      b.min = 0.001, # scaled to 1 but allow overshooting
      b.max = 1.1, # scaled to 1 but allow overshooting
      q.min = rep(0.1,3),  # max value of q , anything larger is not believable
      q.max = rep(2,3),  # max value of q , anything larger is not believable
      bo.mup=rep(-2.5579,3),
      bo.sdp=rep(0.47726,3),
      bp.mup=rep(-2.5579,3),
      bp.sdp=rep(0.47726,3),
      rec.max= c( 10^3, 10^4, 10^2 ),
      K.mu = c( 1.831139,4.170013,0.784308), #for ln
      K.sd = c(0.04595339,0.04350642,0.02571607), #for ln
      r.mu = rep(0.96,3),
      r.sd = rep(0.01041271,3),
      b0.min = c(0.5, 0.5, 0.2),  # prior: mean value possible in  N,S,4X
      b0.max = c(0.8, 0.8, 0.6),  # prior: mean value possible in  N,S,4X
      cv.normal.min = 0.05, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
      cv.normal.max = 0.4, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
      cv.lognormal.min = 0.05, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
      cv.lognormal.max = 0.4, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
  
      IOA = as.matrix(res$B), # observed index of abundance
      IOAcv = as.matrix(res$B.sd ), # observed index of log abundance SD estimates ~ CV
      IREC = as.matrix(res$R), # observed index of abundance
      IRECcv = as.matrix(res$R.sd ), # observed index of log abundance SD estimates ~CV
      CAT = as.matrix(res$L) , # catches  , assume 20% handling mortality and illegal landings
      CAT.min = apply( res$L, 2, min, na.rm=T),
      CAT.max = apply( res$L, 2, max, na.rm=T),
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
    sb$IRECcv[ cfa.north.bad.data, cfa.north ] = mean( sb$IRECcv[,cfa.north ] ,na.rm=T )
    sb$IREC[ cfa.north.bad.data, cfa.north ] = mean( sb$IREC[,cfa.north ] ,na.rm=T )
    cfa.south =  2 # column index
    cfa.south.bad.data = which( as.numeric(rownames(sb$IOA)) <= 1998 )
    sb$IOA[ cfa.south.bad.data, cfa.south ] = NA
    sb$IOAcv[ cfa.south.bad.data, cfa.south ] = mean( sb$IOAcv[,cfa.south ] ,na.rm=T )
    cfa.nodata =   which( as.numeric(rownames(sb$IOA)) <= 2003 )
    sb$IOA[ cfa.nodata , sb$cfa4x ] = NA
    sb$IOAcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IOAcv[,sb$cfa4x] ,na.rm=T )

    # recruitment index has no process model right now ... use min value
    sb$IREC[ cfa.nodata , sb$cfa4x ] = min( sb$IREC[,sb$cfa4x] ,na.rm=T )
    sb$IRECcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IRECcv[,sb$cfa4x] ,na.rm=T )
    sb$IREC = log( sb$IREC )
    
    sb$tomonitor = c( "r", "K", "q", "qs", "r.mu", "r.sd", "b","bp.sd", "bo.sd", "b0", "b0.sd", "F", "TAC",  "C", "P", "B", "rem", "rem.sd", "rem.mu", "REM", "MSY", "BMSY", "FMSY", "Fcrash", "Bdrop", "BX2MSY" )
    sb$jagsmodelname = "biomassdynamic_nonhyper_2014.bugs"

    return(sb)
  }

  if ( DS %in% c("jags", "jags.2015")  ) {
    #AMC: hard coded landings for sens 2004-2014
    res$L[9:19,2] <- c(8.022, 6.407, 4.486,4.942,8.253,10.645,13.150,12.135,11.733,11.309,11.267)

    #AMC: hard coded landings for cfa4x 2004-2014
    res$L[1:9,3] <- c(.004,res$L[1:8,3])

    ###  all data follow this sequence: c("cfanorth", "cfasouth", "cfa4x")
    sb = list(
      b.min = 0.001, # scaled to 1 but allow overshooting
      b.max = 1.1, # scaled to 1 but allow overshooting
      q.min = rep(0.001,3),
      q.max = rep(1,3),
      bo.mup=rep(-2.5579,3),
      bo.sdp=rep(0.47726,3),
      bp.mup=rep(-2.5579,3),
      bp.sdp=rep(0.47726,3),
      rec.max= c( 10^3, 10^4, 10^2 ),
      K.mu = c( 1.831139,4.170013,0.784308), #for ln
      K.sd = c(0.04595339,0.04350642,0.02571607), #for ln
      r.mu = rep(0.96,3),
      r.sd = rep(0.01041271,3),
      b0.min = c(0.5, 0.5, 0.2),  # prior: mean value possible in  N,S,4X
      b0.max = c(0.8, 0.8, 0.6),  # prior: mean value possible in  N,S,4X
      cv.normal.min = 0.05, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
      cv.normal.max = 0.4, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
      cv.lognormal.min = 0.05, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
      cv.lognormal.max = 0.4, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
    # for lognormal: cv = sqrt(exp(sigma^2) - 1); or sigma = sqrt(log(cv^2+ 1) ) ==> sigma = sqrt( log(0.25^2 + 1)) = 0.246 ~ cv -- i.e. cv ~ sd
      IOA = as.matrix(res$B), # observed index of abundance
      IOAcv = as.matrix(res$B.sd ), # observed index of log abundance SD estimates ~ CV
      IREC = as.matrix(res$R), # observed index of abundance
      IRECcv = as.matrix(res$R.sd ), # observed index of log abundance SD estimates ~CV
      CAT = as.matrix(res$L) , # catches  , assume 20% handling mortality and illegal landings
      CAT.min = apply( res$L, 2, min, na.rm=T),
      CAT.max = apply( res$L, 2, max, na.rm=T),
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
    sb$IRECcv[ cfa.north.bad.data, cfa.north ] = mean( sb$IRECcv[,cfa.north ] ,na.rm=T )
    sb$IREC[ cfa.north.bad.data, cfa.north ] = mean( sb$IREC[,cfa.north ] ,na.rm=T )

    cfa.south =  2 # column index
    cfa.south.bad.data = which( as.numeric(rownames(sb$IOA)) <= 1998 )
    sb$IOA[ cfa.south.bad.data, cfa.south ] = NA
    sb$IOAcv[ cfa.south.bad.data, cfa.south ] = mean( sb$IOAcv[,cfa.south ] ,na.rm=T )

    cfa.nodata =   which( as.numeric(rownames(sb$IOA)) <= 2003 )
    sb$IOA[ cfa.nodata , sb$cfa4x ] = NA
    sb$IOAcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IOAcv[,sb$cfa4x] ,na.rm=T )

    # recruitment index has no process model right now ... use min value
    sb$IREC[ cfa.nodata , sb$cfa4x ] = min( sb$IREC[,sb$cfa4x] ,na.rm=T )
    sb$IRECcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IRECcv[,sb$cfa4x] ,na.rm=T )

    sb$IREC = log( sb$IREC )

    sb$tomonitor = c( "r", "K", "q", "qs", "r.mu", "r.sd", "K.sd", "b", "bp.sd", "bo.sd", "b0", "b0.sd", "rem", "rem.sd", "rem.mu", "REM", "MSY", "BMSY", "FMSY", "Fcrash", "Bdrop", "BX2MSY", "F", "TAC",  "C", "P", "B" )
    sb$jagsmodelname = "biomassdynamic_nonhyper_2014.bugs"

    return(sb)

  }
  

 # ---------------------------------

  if (DS %in% c("jags.2016")) {
    sb = list(
      b.min = 0.001, # scaled to 1 but allow overshooting
      b.max = 1.1, # scaled to 1 but allow overshooting
      q.min = 0.001,
      q.max = 2,
      q.mu = rep(1,3),
      q.sd = rep(0.3,3),
      bo.mup=rep(-2.5579,3),
      bo.sdp=rep(0.47726,3),
      bo.min=rep(0,3),
      bo.max=rep(5,3),
      rec.max= c( 10^3, 10^4, 10^2 ),
      bp.mup=rep(-2.5579,3),
      bp.sdp=rep(0.47726,3),
      bp.min=rep(0,3),
      bp.max=rep(5,3),
      rec.max= c( 10^3, 10^4, 10^2 ),
      K.mu = c( 1.831139,4.170013,0.784308), #for ln
      K.sd = c(0.04595339,0.04350642,0.02571607), #for ln
      r.mu = rep(0.96,3),
      r.sd = rep(0.01041271,3),
      b0.min = c(0.5, 0.5, 0.2),  # prior: mean value possible in  N,S,4X
      b0.max = c(0.8, 0.8, 0.6),  # prior: mean value possible in  N,S,4X
      cv.normal.min = 0.05, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
      cv.normal.max = 0.4, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
      cv.lognormal.min = 0.05, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
      cv.lognormal.max = 0.4, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
    # for lognormal: cv = sqrt(exp(sigma^2) - 1); or sigma = sqrt(log(cv^2+ 1) ) ==> sigma = sqrt( log(0.25^2 + 1)) = 0.246 ~ cv -- i.e. cv ~ sd
      IOA = as.matrix(res$B), # observed index of abundance
      IOAcv = as.matrix(res$B.sd ), # observed index of log abundance SD estimates ~ CV
      #IREC = as.matrix(res$R), # observed index of abundance
      #IRECcv = as.matrix(res$R.sd ), # observed index of log abundance SD estimates ~CV
      CAT = as.matrix(res$L) , # catches  , assume 20% handling mortality and illegal landings
      CAT.min = apply( res$L, 2, min, na.rm=T),
      CAT.max = apply( res$L, 2, max, na.rm=T),
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

    sb$tomonitor = c( "r", "K", "q", "qs", "r.mu", "r.sd", "b","bp.sd", "bo.sd","b0", "b0.sd", "rem", "rem.sd", "rem.mu","REM", "MSY", "BMSY", "FMSY", "Fcrash", "Bdrop", "BX2MSY", "F", "TAC",  "C", "P", "B" )

    sb$jagsmodelname = "biomassdynamic_nonhyper_2016.bugs"
    return(sb)
  }


  if (DS %in% c("LaplacesDemon.debug")) {

    # single region test
    K0est =list()
    K0est[["cfanorth"]] = 5
    K0est[["cfasouth"]] = 65
    K0est[["cfa4x"]] = 5
        
    yrs = as.numeric(rownames( res$B))
    nforecasts = 5
    
    sb = list(
      Ndata = length(yrs),
      Nforecasts = nforecasts, # no years for projections
      O = res$B[[debug.region]], # observed index of abundance
      log_O0 = mean( log(res$B[[debug.region]]), na.rm=TRUE ),
      Omax = max( res$B[[debug.region]] *1.5 , na.rm=TRUE),
      removals = res$L[[debug.region]] , # removalsches  , assume 20% handling mortality and illegal landings
      log_R0 = log(mean(res$L[[debug.region]], na.rm=TRUE )/K0est[[debug.region]]),
      er = 0.2,  # target exploitation rate
      ty = which(yrs==2004) ,  # index of the transition year (2004) between spring and fall surveys
      log_r0= log(1),
      log_K0= log(K0est[[debug.region]]),
      log_q0= log(mean(res$B[[debug.region]], na.rm=TRUE)/K0est[[debug.region]]),
      S0= 0.6, # normalised
      cv = 0.4,
      smax =1.25,
      eps = 1e-6, 
      eps_1 = 1 - 1e-6
    )

    return(sb)
  }


}
