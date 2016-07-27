
surplus.production.simple.laplacesdemon.setup = function(Data) {
  
  # zeta = as.vector(mvtnorm::rmvnorm(1, rep(0,Data$S), Sigma, method="chol" ))

  # set up model for a simple surplus production model 
  # to be solved by the Rlibrary: LaplacesDemon (or alternately via penalized Maximum Likelihood) 

  if (0) {
    # to install packages /debug
    install_github("LaplacesDemonR/LaplacesDemon")
    install_github("LaplacesDemonR/LaplacesDemonCpp")
    require(LaplacesDemon)
    require(LaplacesDemonCpp)
  }

  if (is.null( Data$Missing )) {
    mO = which( !is.finite(Data$O))
    mr = which( !is.finite(Data$removals) )
    Data$Missing = list(
      O = mO,
      nO = length(mO) ,
      removals = mr , 
      nremovals = length( mr ),
      n = length( mO) + length(mr) 
    )
  }

  if (is.null( Data$PGF )) {
    # Parameter Generating Function
    Data$PGF = function(Data) {
      r=rnorm(1, log(Data$r0), sd=Data$cv )
      K=rnorm(1, log(Data$K0), sd=Data$cv ) # log scale
      q=rnorm(1, log(Data$q0), sd=Data$cv )
      S_sd=runif( 1 )
      O_sd=runif( 1 )
      S=rnorm( Data$N, log(Data$S0), sd=Data$cv ) # log scale 
      S0=rnorm(1, log(Data$S0), sd=Data$cv ) # log scale
      out = c(r, K, q, S_sd, O_sd, S, S0)
      if (Data$Missing$nO > 0) {
        Omissing = rnorm( Data$Missing$nO, log(Data$Omissing0), sd=Data$cv )
        out = c( out, Omissing )
      }
      if (Data$Missing$nremovals > 0) {
        removalsmissing = rnorm( Data$Missing$nremovals, log(Data$removalsmissing0), sd=Data$cv )
        out = c( out, removalsmissing )
      }
      return( out  )
    }
    Data$PGF = compiler::cmpfun(Data$PGF)
  }

  if (is.null( Data$parm.names )) {
    # paramater names and initial values
    pnames = list(
      r = Data$r0,
      K = Data$K0,
      q = Data$q0,
      S_sd = Data$cv,
      O_sd = Data$cv,
      S = rep( Data$S0, Data$N),
      S0 = Data$S0
    )
    if (Data$Missing$nO > 0) pnames$Omissing = rep( Data$Omissing0, Data$Missing$nO )
    if (Data$Missing$nremovals > 0)  pnames$removalsmissing = rep( Data$removalsmissing0, Data$Missing$nremovals )
    Data$parm.names = as.parm.names( pnames )
  }

  if (is.null( Data$idx )) {
    # index position of paramaters
    Data$idx = list(
      r=grep("\\<r\\>", Data$parm.names),
      K=grep("\\<K\\>", Data$parm.names),
      q=grep("\\<q\\>", Data$parm.names),
      S_sd=grep("\\<S_sd\\>", Data$parm.names),
      O_sd=grep("\\<O_sd\\>", Data$parm.names),
      S=grep("\\<S\\>", Data$parm.names),
      S0=grep("\\<S0\\>", Data$parm.names)
    )
    if (Data$Missing$nO > 0) Data$idx$Omissing = grep("\\<Omissing\\>", Data$parm.names)
    if (Data$Missing$nremovals > 0)  Data$idx$removalsmissing = grep("\\<removalsmissing\\>", Data$parm.names)
  }


  if (is.null( Data$mon.names )) {
    Data$mon.names = c("LP", "r", "K", "q" )
    # Data$mon.names = c("LP", "r", "K", "q", paste0("S",1:Data$N), paste0("AR",1:(Data$N-1) ) )
  }


  if (is.null( Data$Model )) {
    Data$Model = function(parm, Data) {
      
      Spred = Opred = rep(0, Data$N )   # initialize a few storage vectors

      parm_exp = exp(parm) # all parms in log-space to keep positive valued .. return to correct scale
      # NOTE: for lognormal: cv = sqrt(exp(sigma^2) - 1); 
      # or sigma = sqrt(log(cv^2+ 1) ) ==> sigma = sqrt( log(0.25^2 + 1)) = 0.246 ~ cv -- i.e. cv ~ sd
      q = parm_exp[Data$idx$q]
      r = parm_exp[Data$idx$r]
      K =  parm_exp[Data$idx$K] 
      Spred[1] =  parm_exp[Data$idx$S0] 
      S =  parm_exp[Data$idx$S] 
      O_sd = parm_exp[Data$idx$O_sd] 
      S_sd = parm_exp[Data$idx$S_sd] 
      
      loglik = c()
      loglik[Data$parm.names] = 0
      loglik[Data$idx$q] = dnorm( parm[Data$idx$r], log(Data$q0), Data$cv, TRUE ) ;
      loglik[Data$idx$r] = dnorm( parm[Data$idx$r], log(Data$r0), Data$cv, TRUE ) ;
      loglik[Data$idx$K] = dnorm( parm[Data$idx$K], log(Data$K0), Data$cv, TRUE ) ;
      loglik[Data$idx$S0] = dnorm( parm[Data$idx$S0], log(Data$S0), Data$cv, TRUE ) ;
      loglik[Data$idx$O_sd] = LaplacesDemonCpp::dhalfcauchy( parm[Data$idx$O_sd], Data$cv, TRUE );
      loglik[Data$idx$S_sd] = LaplacesDemonCpp::dhalfcauchy( parm[Data$idx$S_sd], Data$cv, TRUE );

      R = Data$removals/K ;
      if (exists( "Missing", Data ) ) {
        if (Data$Missing$nremovals > 0) {
          R[Data$Missing$removals ] = rlnorm( Data$Missing$nremovals, 
            parm[Data$idx$removalsmissing], sd=Data$cv )
          R = .Internal(pmax(na.rm=TRUE, R, Data$eps )) ;
          loglik[Data$idx$removalsmissing] = dnorm( parm[Data$idx$removalsmissing], 
            log(R[Data$Missing$removals ]), Data$cv )
        }
      }

      AR = 1.0 + r*(1-S[Data$jj]) ;  # "autocorelation" component
      Spred[ Data$ii ] = S[Data$jj] * AR - R[Data$jj] ;  # simple logistic
      Spred = bio.utilities::truncate.vector( Spred, lower=Data$eps ) ;
      loglik[Data$idx$S] = dnorm( parm[Data$idx$S], log(Spred), S_sd, TRUE );

      # Likelihoods for observation model
      i = 1 ;                   Opred[i] = S[i] - R[i] ;
      i = 2:(Data$ty-1);        Opred[i] = S[i] - R[i-1] ;
      i = Data$ty;              Opred[i] = S[i] - (R[i-1] + R[i])/2 ;
      i = (Data$ty+1):Data$N ;  Opred[i] = S[i] - R[i] ;
      Opred = bio.utilities::truncate.vector( K*q*Opred, lower=Data$eps )  ;
      if (exists( "Missing", Data ) ) {
        if (Data$Missing$nO > 0) {
          Data$O[Data$Missing$O ] = rlnorm( Data$Missing$nO, parm[Data$idx$Omissing], sd=Data$cv )
          Data$O[Data$Missing$O ] = .Internal(pmax(na.rm=TRUE, Data$O[Data$Missing$O ], Data$eps )) ;
          loglik[Data$idx$Omissing] = dnorm( parm[Data$idx$Omissing], log(Data$O[Data$Missing$O ]), Data$cv )
        }
      }
      ll_obs = dnorm( log(Data$O), log(Opred), O_sd, TRUE )
      
      Smon = Rmon = ER = F = B = C = rep(0, Data$MN )
      Smon[1:Data$N] = S
      Rmon[1:Data$N] = R
      for( i in (Data$N+1):Data$MN ){
        Smon[i] = Smon[i-1] * (1.0 + r*(1.0-Smon[i-1]) ) - Rmon[i-1]
        Rmon[i] = Smon[i-1] * Data$er 
      }
      Smon = bio.utilities::truncate.vector( Smon, lower=Data$eps ) 
      Rmon = bio.utilities::truncate.vector( Rmon, lower=Data$eps ) 
      
      # monitoring
      ER = Rmon / Smon ;
      ER = bio.utilities::truncate.vector( ER, Data$eps, 1-Data$eps ) 
      
      B = Smon*K
      C = Rmon*K
      F = -log( 1 - ER) ; # fishing mortality
            
      LL = sum( ll_obs )  # log likelihood (of the data)
      LP = LL + sum( loglik )  # log posterior
      out = list( LP=LP, Dev=-2*LL, Monitor=c(LP, r, K, q), yhat=S*K, parm=parm )
      return( out  )
    }
    Data$Model.ML  = compiler::cmpfun( function(...) (Data$Model(...)$Dev / 2) )  # i.e. - log likelihood
    Data$Model.PML = compiler::cmpfun( function(...) (- Data$Model(...)$LP) ) #i.e., - log posterior 
    Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed
  }
  return (Data)

}

