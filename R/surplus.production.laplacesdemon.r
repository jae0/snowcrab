
surplus.production.laplacesdemon = function(Data) {

# set up model for a simple surplus production model 
# to be solved by the Rlibrary: LaplacesDemon (or alternately via Maximum Likelihood) 

  if (0) {
    # to install packages /debug
    install_github("lazycrazyowl/LaplacesDemon")
    install_github("lazycrazyowl/LaplacesDemonCpp")
    require(LaplacesDemon)
    require(LaplacesDemonCpp)
  }

  if (is.null( Data$PGF )) {
    # Parameter Generating Function
    Data$PGF = function(Data) {
      r=rnorm(1, Data$r0, sd=Data$cv )
      K=rnorm(1, Data$K0, sd=Data$cv ) # log scale
      q=rnorm(1, Data$q0, sd=Data$cv )
      S_sd=runif( Data$N )
      O_sd=runif( Data$N )
      S=rnorm( Data$N, Data$S0, sd=Data$cv )
      S0=rnorm(1, Data$S0, sd=Data$cv )
      return( c(r, K, q, S_sd, O_sd, S, S0) )
    }
  }

  if (is.null( Data$parm.names )) {
    # paramater names and initial values
    Data$parm.names = as.parm.names( list(
      r = Data$r0,
      K = Data$K0,
      q = Data$q0,
      S_sd = rep( Data$cv, Data$N),
      O_sd = rep( Data$cv, Data$N),
      S = rep( Data$S0, Data$N),
      S0 = Data$S0
    ))
  }

  if (is.null( Data$Data$pos )) {
    # index position of paramaters
    Data$pos = list(
      r=grep("\\<r\\>", Data$parm.names),
      K=grep("\\<K\\>", Data$parm.names),
      q=grep("\\<q\\>", Data$parm.names),
      S_sd=grep("\\<S_sd\\>", Data$parm.names),
      O_sd=grep("\\<O_sd\\>", Data$parm.names),
      S=grep("\\<S\\>", Data$parm.names),
      S0=grep("\\<S0\\>", Data$parm.names)
    )
  }

  if (is.null( Data$Data$pos )) {
    Data$mon.names = c("LP", "r", "K", "q", paste0("S",1:Data$N), paste0("AR",1:(Data$N-1) ) )
  }

  if (is.null( Data$Data$pos )) {
    Data$Model = function(parm, Data) {
      Spred = Opred = ER = F = rep(0, Data$N )   # initialize
      # Priors .. remember that these are on log-scale
      loglik = c()
      loglik[Data$parm.names] = 0
      loglik[Data$pos$q] = dnorm( parm[Data$pos$q], Data$q0, Data$cv,  TRUE ) ;
      loglik[Data$pos$r] = dnorm( parm[Data$pos$r], Data$r0, Data$cv, TRUE ) ;
      loglik[Data$pos$K] = dnorm( parm[Data$pos$K], Data$K0, Data$cv, TRUE ) ;
      loglik[Data$pos$S0] = dnorm( parm[Data$pos$S0], Data$S0, Data$cv, TRUE ) ;

      # Parameters: back transform to force positive only solutions
      r = exp( parm[Data$pos$r] );
      K = exp( parm[Data$pos$K] );
      q = exp( parm[Data$pos$q] );
      Spred[1] = exp( parm[Data$pos$S0] );
      S = exp( parm[Data$pos$S] ) ;
      O_sd = exp(parm[Data$pos$O_sd])
      S_sd = exp(parm[Data$pos$S_sd])

      # these are SD's on log scale (0,1) is sufficient
      # continue with SD priors as these are now on correct scale
      loglik[Data$pos$O_sd] = dnorm( O_sd, Data$cv, 0.1, TRUE );
      loglik[Data$pos$S_sd] = dnorm( S_sd, Data$cv, 0.1, TRUE );

      # Likelihoods for process model
      R = Data$removals/K ;
      AR = 1.0 + r*(1-S[Data$jj]) ;  # autocorelation component
      Spred[ Data$ii ] = S[Data$jj] * AR - R[Data$jj] ;  # simple logistic
      Spred = .Internal(pmax(na.rm=TRUE, Spred, Data$eps )) ;
      loglik[Data$pos$S] = dnorm( log(S), log(Spred), S_sd, TRUE );

      # Likelihoods for observation model
      i = 1 ;                   Opred[i] = S[i] - R[i] ;
      i = 2:(Data$ty-1);        Opred[i] = S[i] - R[i-1] ;
      i = Data$ty;              Opred[i] = S[i] - (R[i-1] + R[i])/2 ;
      i = (Data$ty+1):Data$N ;  Opred[i] = S[i] - R[i] ;
      Opred = .Internal(pmax( na.rm=TRUE, K*q*Opred, Data$eps ) ) ;
      ll_obs = dnorm( log(O), log(Opred), O_sd, TRUE )

      # Removals
    #  for( int i=1; i<nt; i++){
    #    nloglik[2] -= dnorm( R(i-1), Rp(i-1), Rp(i-1)*cv_R, true );
    #  }

    # forecast and monitoring
    #  for( i 1:Data$N ){
    #    ER[i] = R[i] / S[i] ;
      #  B(i) <- S(i)*K
      #  C(i) <- R(i)*K
    #    F[i] = -log(1 - ER[i] ) ; # fishing mortality
    #  }

      LL = sum( ll_obs )  # log likelihood (of the data)
      LP = LL + sum( loglik )  # log posterior
      out = list( LP=LP, Dev=-2*LL, Monitor=c(LP, r, K, q, AR), yhat=S*K, parm=parm )
      return( out  )
    }
    Data$Model <- compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed
  }
  return (Data)

}
