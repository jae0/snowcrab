


surplus.production.simple.laplacesdemon.setup = function(Data) {
  
  # set up model for a simple surplus production model 
  # to be solved by the Rlibrary: LaplacesDemon (or alternately via penalized Maximum Likelihood) 


  # ----------------------------------
  # Identify location and number of missing values -- prediction locations are treated the same way

  # compute a few things here that are constant in the model
  Data$N = Data$Ndata + Data$Nforecasts # no years with data + projections
  Data$q0 = exp( Data$log_q0) 
  Data$r0 = exp( Data$log_r0) 
  Data$K0 = exp( Data$log_K0) 
  Data$removals0 = exp(Data$log_removals0)
  Data$O0 = exp(Data$log_O0)
  Data$O_range = range( Data$O, na.rm=TRUE)

  # ----------------------------------

  Data$Missing = list(
    O = list( 
      n=length(intersect( 1:Data$Ndata, which( !is.finite(Data$O) ) )), 
      idx=intersect( 1:Data$Ndata, which( !is.finite(Data$O) ) ) 
    ), # observations of survey index
    removals = list( 
      n=length(intersect( 1:Data$Ndata, which( !is.finite(Data$removals) ) ) ), 
      idx=intersect( 1:Data$Ndata, which( !is.finite(Data$removals) ) )  
    )  # observations of removals
  )

  # ----------------------------------

  if ( Data$Nforecasts > 0 ) {
    # extend data for forecasts
    Data$O = c( Data$O, rep(NA, Data$Nforecasts) )
    Data$removals = c( Data$removals, rep(NA, Data$Nforecasts) )
    Data$Forecasts = list( 
      n = Data$Nforecasts,
      idx = Data$Ndata + 1:Data$Nforecasts,
      idx_1 = (Data$Ndata + 1:Data$Nforecasts) -1 
    )
  }

  # ----------------------------------
  # misc indexing
  Data$idx = list(
    t1 = 1,
    t2 = 2:(Data$ty-1),
    t2_1 = (2:(Data$ty-1) ) - 1,
    t3 =  Data$ty,
    t3_1 =  Data$ty-1,
    t4 =  (Data$ty+1):Data$N,
    icurrent  = 2:Data$N,
    iprevious = 1:(Data$N-1),
    idata = 1:Data$Ndata
  )

  # ----------------------------------
  # paramater names and initial values
  Data$parm.names = as.parm.names( list(
    r = Data$r0,
    K = Data$K0,
    q = Data$q0,
    r_sd = Data$r0,
    K_sd = Data$K0,
    q_sd = Data$q0,
    S0_sd = Data$S0,
    S_sd = Data$S0,
    O_sd = Data$O0,
    R_sd = Data$removals0,
    S = rep( Data$S0, Data$N),
    S0 = Data$S0
  ) )

  # ----------------------------------
  # index position of paramaters
  Data$pos = list(
    r=grep("\\<r\\>", Data$parm.names),
    K=grep("\\<K\\>", Data$parm.names),
    q=grep("\\<q\\>", Data$parm.names),
    r_sd = grep("\\<r_sd\\>", Data$parm.names), 
    K_sd = grep("\\<K_sd\\>", Data$parm.names),
    q_sd = grep("\\<q_sd\\>", Data$parm.names),
    S0_sd = grep("\\<S0_sd\\>", Data$parm.names),
    S_sd=grep("\\<S_sd\\>", Data$parm.names),
    O_sd=grep("\\<O_sd\\>", Data$parm.names),
    R_sd =grep("\\<R_sd\\>", Data$parm.names),
    S=grep("\\<S\\>", Data$parm.names),
    S0=grep("\\<S0\\>", Data$parm.names)
  )
 
  # ----------------------------------
  # monitoring nodes
  Data$mon.names = c("LP", "r", "K", "q" )
  # Data$mon.names = c("LP", "r", "K", "q", paste0("S",1:Data$N), paste0("AR",1:(Data$N-1) ) )


  # ----------------------------------
  # Parameter Generating Function
  # these parameters are operate on the log-scale to force positive values ...
  Data$PGF = function(Data) {
    r = runif( 1, Data$r0/2, Data$r0*2 )
    K = runif( 1, Data$K0/2, Data$K0*2 ) 
    q = runif( 1, Data$q0/2, Data$q0*2 )
    r_sd = runif( 1, Data$eps, Data$r0/4)
    K_sd = runif( 1, Data$eps, Data$K0/4)
    q_sd = runif( 1, Data$eps, Data$q0/4)
    S0_sd = runif( 1, Data$eps, Data$S0/4)
    S_sd = runif( 1, Data$eps, 0.5 );
    O_sd = runif( 1, Data$eps, Data$O0/4 );
    R_sd = runif( 1, Data$eps, Data$removals0/4 );
    S =  rbeta( Data$N, 5, 5 ) + Data$eps  
    S0 = runif( 1, 0.1, 0.8 ); 
    out = log(c(r, K, q, r_sd, K_sd, q_sd, S0_sd, S_sd, O_sd, R_sd, S, S0) )
    return( out  )
  }
  Data$PGF = compiler::cmpfun(Data$PGF)


  # ----------------------------------
  # define the model that generates the loglikelihoods (LL) and the log-posteriors (LP)
  Data$Model = function(parm, Data) {
    
    Spred = Opred = rep(0, Data$N )   # initialize a few storage vectors
      # NOTE: for lognormal: cv = sqrt(exp(sigma^2) - 1); 
    # or sigma = sqrt(log(cv^2+ 1) ) ==> sigma = sqrt( log(0.25^2 + 1)) = 0.246 ~ cv -- i.e. cv ~ sd
    pm = exp(parm)

    q = pm[Data$pos$q]
    r = pm[Data$pos$r]
    K = pm[Data$pos$K]
    S = pm[Data$pos$S] 
    Spred[1] = pm[Data$pos$S0] 
   
   # SD params stay on log-scale
    q_sd = pm[Data$pos$q_sd]
    r_sd = pm[Data$pos$r_sd]
    K_sd = pm[Data$pos$K_sd]
    S0_sd = pm[Data$pos$S0_sd]
    O_sd = pm[Data$pos$O_sd]  
    S_sd = pm[Data$pos$S_sd] 
    R_sd = pm[Data$pos$r_sd] 
    
    llkimpute = 0
    llkprior = c()
    llkprior[Data$parm.names] = 0
    llkprior[Data$pos$q] = dnorm( parm[Data$pos$q], Data$log_q0, q_sd, TRUE ) ;
    llkprior[Data$pos$r] = dnorm( parm[Data$pos$r], Data$log_r0, r_sd, TRUE ) ;
    llkprior[Data$pos$K] = dnorm( parm[Data$pos$K], Data$log_K0, K_sd, TRUE ) ;
    llkprior[Data$pos$S0] = dnorm( parm[Data$pos$S0], log(Data$S0), S0_sd, TRUE ) ;

    llkprior[Data$pos$q_sd] = dhalfcauchy( pm[Data$pos$q_sd], Data$q0, TRUE );
    llkprior[Data$pos$r_sd] = dhalfcauchy( pm[Data$pos$r_sd], Data$r0, TRUE );
    llkprior[Data$pos$K_sd] = dhalfcauchy( pm[Data$pos$K_sd], Data$K0, TRUE );
    llkprior[Data$pos$S0_sd] = dhalfcauchy( pm[Data$pos$S0_sd], Data$S0, TRUE );
    llkprior[Data$pos$O_sd] = dhalfcauchy( pm[Data$pos$O_sd], Data$O0, TRUE );
    llkprior[Data$pos$S_sd] = dhalfcauchy( pm[Data$pos$S_sd], 0.5, TRUE );
    llkprior[Data$pos$R_sd] = dhalfcauchy( pm[Data$pos$R_sd], Data$removals0, TRUE );

    R = Data$removals/K ;# make sure it is producing sensible values:
    # impute missing data 

    if ( Data$Missing$removals$n > 0 ) {
      R[Data$Missing$removals$idx ] = runif( Data$Missing$removals$n, Data$eps, Data$K0 )
      llkimpute = llkimpute + sum( dnorm( log(R[Data$Missing$removals$idx ]), Data$log_removals0, R_sd, TRUE ) );
    }
    if ( Data$Missing$O$n > 0 ) {
      Data$O[Data$Missing$O$idx ] = runif( Data$Missing$O$n, Data$eps, Data$K0 )
      llkimpute = llkimpute + sum(dnorm( log(Data$O[Data$Missing$O$idx ]), Data$log_O0, O_sd, TRUE ) ) ;    
    }
    if ( Data$Forecasts$n > 0 ) {
      R[Data$Forecasts$idx] = S[Data$Forecasts$idx_1] * Data$er 
      llkimpute = llkimpute + sum(dnorm( log(R[Data$Forecasts$idx]), Data$log_removals0, R_sd, TRUE )) ;
      Data$O[Data$Forecasts$idx] = runif( Data$Forecasts$n, Data$O_range[1], Data$O_range[2]  )
      # likelihood below
    }

    # process model -- simple logistic
    Spred[ Data$idx$icurrent] = S[Data$idx$iprevious] * 1.0 + r*{1-S[Data$idx$iprevious]} - R[Data$idx$iprevious] ;  # 
    Spred = bio.utilities::truncate.vector( Spred, lower=Data$eps ) 
    llkprior[Data$pos$S] = dnorm( parm[Data$pos$S] , log(Spred), S_sd, TRUE ) 

    # observation model
    Opred[Data$idx$t1] = S[Data$idx$t1] - R[Data$idx$t1] ; # first year approximation
    Opred[Data$idx$t2] = S[Data$idx$t2] - R[Data$idx$t2_1] ;
    Opred[Data$idx$t3] = S[Data$idx$t3] - (R[Data$idx$t3_1] + R[Data$idx$t3])/2 ; # transition year from Spring to Autumn survey
    Opred[Data$idx$t4] = S[Data$idx$t4] - R[Data$idx$t4] ;
    Opred = K*q*Opred
    Opred = bio.utilities::truncate.vector( Opred, Data$eps )
    llkdata = dnorm( log(Data$O), log(Opred), O_sd, TRUE )
 
    # additional computed variables of interest 
    ER = R / S ;
    ER = bio.utilities::truncate.vector( ER, Data$eps, Data$eps_1 ) 
    B = S*K
    C = R*K
    F = -log( 1 - ER) ; # fishing mortality
    AR = S[ Data$idx$icurrent] / S[Data$idx$iprevious] ;  # "autocorelation" component

    LL = sum( llkdata )  # log likelihood (of the data)
    LP = LL + sum( llkprior ) + sum(llkimpute) # log posterior
    out = list( LP=LP, Dev=-2*LL, Monitor=c(LP, r, K, q), yhat=S*K, parm=parm )
    return( out  )
  }

  Data$Model.ML  = compiler::cmpfun( function(...) (Data$Model(...)$Dev / 2) )  # i.e. - log likelihood
  Data$Model.PML = compiler::cmpfun( function(...) (- Data$Model(...)$LP) ) #i.e., - log posterior 
  Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed

  return (Data)

}

