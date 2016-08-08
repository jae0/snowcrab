


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
  Data$R0 = exp(Data$log_R0)
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
  Data$parm.names = list(
    r = Data$r0,
    K = Data$K0,
    q = Data$q0,
    r_sd = Data$r0,
    K_sd = Data$K0,
    q_sd = Data$q0,
    S0_sd = Data$S0,
    S_sd = Data$S0,
    O_sd = Data$O0,
    R_sd = Data$R0,
    S = rep( Data$S0, Data$N),
    S0 = Data$S0
  ) 
  if (Data$Missing$removals$n > 0) Data$parm.names$removals_imp = rep(0, Data$Missing$removals$n)
  if (Data$Missing$O$n > 0) Data$parm.names$O_imp = rep(0, Data$Missing$O$n)
  Data$parm.names = as.parm.names( Data$parm.names )

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
  if (Data$Missing$removals$n > 0) Data$pos$removals_imp = grep("\\<removals_imp\\>", Data$parm.names)
  if (Data$Missing$O$n > 0) Data$pos$O_imp = grep("\\<O_imp\\>", Data$parm.names)

  # ----------------------------------
  # monitoring nodes
  Data$mon.names = c("LP", "r", "K", "q" )
  # Data$mon.names = c("LP", "r", "K", "q", paste0("S",1:Data$N), paste0("AR",1:(Data$N-1) ) )


  # ----------------------------------
  # Parameter Generating Function
  # these parameters are operate on the log-scale to force positive values ...
  Data$param$means = c( 
    r=Data$r0, K=Data$K0, q=Data$q0, 
    r_sd=Data$r0/2, K_sd=Data$K0/2,  q_sd=Data$q0/2, 
    S0_sd=Data$S0/2, S_sd=0.2, O_sd=Data$O0/2, R_sd=Data$R0/2,
    S=rep( 0.5, Data$N), S0=Data$S0
  )
  if (Data$Missing$removals$n > 0) {
    Data$param$means = c( Data$param$means, removals_imp=rep(Data$R0, Data$Missing$removals$n ) )
  }
  if (Data$Missing$O$n > 0) {
    Data$param$means = c( Data$param$means, O_imp=rep(Data$O0, Data$Missing$O$n) )
  }
  Data$param$sd = Data$param$means / 2 
  Data$param$n = length(Data$param$means)

  Data$PGF = function(Data) {
    out = LaplacesDemonCpp::rnorm(Data$param$n, Data$param$means, Data$param$sd )
    i = which( out < Data$eps )
    if (length(i)>0) out[i] = -out[i] # reflect onto positive side
    return( out  )
  }
  Data$PGF = compiler::cmpfun(Data$PGF)


  # ----------------------------------
  # define the model that generates the loglikelihoods (LL) and the log-posteriors (LP)
  Data$Model = function(parm, Data) {
    
    Spred = Opred = rep(0, Data$N )   # initialize a few storage vectors
      # NOTE: for lognormal: cv = sqrt(exp(sigma^2) - 1); 
    # or sigma = sqrt(log(cv^2+ 1) ) ==> sigma = sqrt( log(0.25^2 + 1)) = 0.246 ~ cv -- i.e. cv ~ sd
    i = which( parm < Data$eps )
    if (length(i)>0) parm[i] = -parm[i] # reflect onto positive side

    q =  parm[Data$pos$q]
    r =  parm[Data$pos$r]
    K =  parm[Data$pos$K]
    S =  parm[Data$pos$S]
    Spred[1] =  parm[Data$pos$S0]
   
   # SD params stay on log-scale
    q_sd = parm[Data$pos$q_sd]
    r_sd = parm[Data$pos$r_sd]
    K_sd = parm[Data$pos$K_sd]
    S0_sd = parm[Data$pos$S0_sd]
    O_sd = parm[Data$pos$O_sd]
    S_sd = parm[Data$pos$S_sd]
    R_sd = parm[Data$pos$R_sd]
    
    llkimpute = 0
    llkprior = c()
    llkprior[Data$parm.names] = 0
    llkprior[Data$pos$q] = dnorm( log(parm[Data$pos$q]), Data$log_q0, q_sd, log=TRUE ) ;
    llkprior[Data$pos$r] = dnorm( log(parm[Data$pos$r]), Data$log_r0, r_sd, log=TRUE ) ;
    llkprior[Data$pos$K] = dnorm( log(parm[Data$pos$K]), Data$log_K0, K_sd, log=TRUE ) ;
    llkprior[Data$pos$S0] = dnorm( log(parm[Data$pos$S0]), log(Data$S0), S0_sd, log=TRUE ) ;
    llkprior[Data$pos$q_sd] = dgamma( q_sd, shape=Data$q0, scale=100, log=TRUE );
    llkprior[Data$pos$r_sd] = dgamma( r_sd, shape=Data$r0, scale=100, log=TRUE );
    llkprior[Data$pos$K_sd] = dgamma( K_sd, shape=Data$K0, scale=100, log=TRUE );
    llkprior[Data$pos$S0_sd] = dgamma( S0_sd, shape=Data$S0, scale=100, log=TRUE );
    llkprior[Data$pos$O_sd] = dgamma( O_sd, shape=Data$O0, scale=100, log=TRUE );
    llkprior[Data$pos$S_sd] = dgamma( S_sd, shape=0.5, scale=100, log=TRUE );
    llkprior[Data$pos$R_sd] = dgamma( R_sd, shape=Data$R0, scale=100, log=TRUE );

    R = Data$removals/K ;# make sure it is producing sensible values:
    # impute missing data 
    if ( Data$Missing$removals$n > 0 ) {
      parm[Data$pos$removals_imp] = R[Data$Missing$removals$idx ] = LaplacesDemonCpp::interval( parm[Data$pos$removals_imp], Data$eps, Data$smax)  
      llkimpute = llkimpute + sum( dnorm( log(R[Data$Missing$removals$idx ]), Data$log_R0, R_sd, log=TRUE ) );
    }

    if ( Data$Missing$O$n > 0 ) {
      parm[Data$pos$O_imp] = Data$O[Data$Missing$O$idx ] = LaplacesDemonCpp::interval( parm[Data$pos$O_imp], Data$eps, Data$K0)  
      llkimpute = llkimpute + sum(dnorm( log(Data$O[Data$Missing$O$idx ]), Data$log_O0, O_sd, log=TRUE ) ) ;    
    }

    if ( Data$Forecasts$n > 0 ) {
      R[Data$Forecasts$idx] = S[Data$Forecasts$idx_1] * Data$er 
      llkimpute = llkimpute + sum(dnorm( log(R[Data$Forecasts$idx]), Data$log_R0, R_sd, log=TRUE )) ;
      Data$O[Data$Forecasts$idx] = runif( Data$Forecasts$n, Data$O_range[1], Data$O_range[2]  )
      # likelihood for Data$O is below
    }

    # process model -- simple logistic
    Spred[ Data$idx$icurrent] = S[Data$idx$iprevious] * 1.0 + r*{1-S[Data$idx$iprevious]} - R[Data$idx$iprevious] ;  # 
    Spred = bio.utilities::truncate.vector( Spred, lower=Data$eps ) 
    llkprior[Data$pos$S] = dnorm( log(parm[Data$pos$S]), log(Spred), S_sd, log=TRUE ) 

    # observation model
    Opred[Data$idx$t1] = S[Data$idx$t1] - R[Data$idx$t1] ; # first year approximation
    Opred[Data$idx$t2] = S[Data$idx$t2] - R[Data$idx$t2_1] ;
    Opred[Data$idx$t3] = S[Data$idx$t3] - (R[Data$idx$t3_1] + R[Data$idx$t3])/2 ; # transition year from Spring to Autumn survey
    Opred[Data$idx$t4] = S[Data$idx$t4] - R[Data$idx$t4] ;
    Opred = K*q*Opred
    Opred = bio.utilities::truncate.vector( Opred, Data$eps )
    llkdata = dnorm( log(Data$O), log(Opred), O_sd, log=TRUE )
 
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

  print (Data$Model( parm=Data$PGF(Data), Data=Data ) ) # test to see if return values are sensible

  return (Data)

}

