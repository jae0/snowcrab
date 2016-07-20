
surplus.production.allregions.laplacesdemon.setup = function(LD) {

  # set up model for a simple surplus production model 
  # to be solved by the Rlibrary: LaplacesDemon (or alternately via penalized Maximum Likelihood) 

  if (0) {
    # to install packages /debug
    install_github("lazycrazyowl/LaplacesDemon")
    install_github("lazycrazyowl/LaplacesDemonCpp")
    require(LaplacesDemon)
    require(LaplacesDemonCpp)
  }

  if (is.null( LD$Missing )) {
    LD$Missing = list(
      O = which( !is.finite(LD$O)),
      nO = length(which( !is.finite(LD$O))) ,
      removals = which( !is.finite(LD$removals) ) , 
      nremovals = length( which( !is.finite(LD$removals) )  )
    )
  }

  if (is.null( LD$PGF )) {
    # Parameter Generating Function
    LD$PGF = function(LD) {
      r=rnorm(1, LD$r0, sd=LD$cv )
      K=rnorm(1, LD$K0, sd=LD$cv ) # log scale
      q=rnorm(1, LD$q0, sd=LD$cv )
      S_sd=runif( LD$N )
      O_sd=runif( LD$N )
      S=rnorm( LD$N, LD$S0, sd=LD$cv )
      S0=rnorm(1, LD$S0, sd=LD$cv )
      out = c(r, K, q, S_sd, O_sd, S, S0)
      if (LD$Missing$nO > 0) {
        Omissing = rnorm( LD$Missing$nO, LD$Omissing0, sd=LD$cv*LD$Omissing0 )
        out = c( out, Omissing )
      }
      if (LD$Missing$nremovals > 0) {
        removalsmissing = rnorm( LD$Missing$nremovals, LD$removalsmissing0, sd=LD$cv*LD$removalsmissing0 )
        out = c( out, removalsmissing )
      }
      return( out  )
    }
    LD$PGF = compiler::cmpfun(LD$PGF)
  }

  if (is.null( LD$parm.names )) {
    # paramater names and initial values
    pnames = list(
      r = LD$r0,
      K = LD$K0,
      q = LD$q0,
      S_sd = rep( LD$cv, LD$N),
      O_sd = rep( LD$cv, LD$N),
      S = rep( LD$S0, LD$N),
      S0 = LD$S0
    )
    if (LD$Missing$nO > 0) pnames$Omissing = rep( LD$Omissing0, LD$Missing$nO )
    if (LD$Missing$nremovals > 0)  pnames$removalsmissing = rep( LD$removalsmissing0, LD$Missing$nremovals )
    LD$parm.names = as.parm.names( pnames )
  }

  if (is.null( LD$i )) {
    # index position of paramaters
    LD$i = list(
      r=grep("\\<r\\>", LD$parm.names),
      K=grep("\\<K\\>", LD$parm.names),
      q=grep("\\<q\\>", LD$parm.names),
      S_sd=grep("\\<S_sd\\>", LD$parm.names),
      O_sd=grep("\\<O_sd\\>", LD$parm.names),
      S=grep("\\<S\\>", LD$parm.names),
      S0=grep("\\<S0\\>", LD$parm.names)
    )
    if (LD$Missing$nO > 0) LD$i$Omissing = grep("\\<Omissing\\>", LD$parm.names)
    if (LD$Missing$nremovals > 0)  LD$i$removalsmissing = grep("\\<removalsmissing\\>", LD$parm.names)
  }


  if (is.null( LD$mon.names )) {
    LD$mon.names = c("LP", "r", "K", "q" )
    # LD$mon.names = c("LP", "r", "K", "q", paste0("S",1:LD$N), paste0("AR",1:(LD$N-1) ) )
  }


  if (is.null( LD$Model )) {
    LD$Model = function(parm, LD) {
      Spred = Opred = ER = F = rep(0, LD$N )   # initialize
      # Priors .. remember that these are on log-scale
      loglik = c()
      loglik[LD$parm.names] = 0
      loglik[LD$i$q] = dnorm( parm[LD$i$q], LD$q0, LD$cv,  TRUE ) ;
      loglik[LD$i$r] = dnorm( parm[LD$i$r], LD$r0, LD$cv, TRUE ) ;
      loglik[LD$i$K] = dnorm( parm[LD$i$K], LD$K0, LD$cv, TRUE ) ;
      loglik[LD$i$S0] = dnorm( parm[LD$i$S0], LD$S0, LD$cv, TRUE ) ;

      # Parameters: back transform to force iitive only solutions
      r = exp( parm[LD$i$r] );
      K = exp( parm[LD$i$K] );
      q = exp( parm[LD$i$q] );
      Spred[1] = exp( parm[LD$i$S0] );
      S = exp( parm[LD$i$S] ) ;
      O_sd = exp(parm[LD$i$O_sd])
      S_sd = exp(parm[LD$i$S_sd])

      # these are SD's on log scale (0,1) is sufficient
      # continue with SD priors as these are now on correct scale
      loglik[LD$i$O_sd] = dnorm( O_sd, LD$cv, 0.1, TRUE );
      loglik[LD$i$S_sd] = dnorm( S_sd, LD$cv, 0.1, TRUE );

      # Likelihoods for process model
      R = LD$removals/K ;
      if (LD$Missing$nremovals > 0) {
        R[LD$Missing$removals ] = rnorm( LD$Missing$nremovals, 
          parm[LD$i$removalsmissing], sd=LD$cv*parm[LD$i$removalsmissing] )
        R = .Internal(pmax(na.rm=TRUE, R, LD$eps )) ;
        loglik[LD$i$removalsmissing] = dnorm( parm[LD$i$removalsmissing], 
          R[LD$Missing$removals ], LD$cv*R[LD$Missing$removals ] )
      }
      AR = 1.0 + r*(1-S[LD$jj]) ;  # autocorelation component
      Spred[ LD$ii ] = S[LD$jj] * AR - R[LD$jj] ;  # simple logistic
      Spred = .Internal(pmax(na.rm=TRUE, Spred, LD$eps )) ;
      loglik[LD$i$S] = dnorm( log(S), log(Spred), S_sd, TRUE );

      # Likelihoods for observation model
      i = 1 ;                   Opred[i] = S[i] - R[i] ;
      i = 2:(LD$ty-1);        Opred[i] = S[i] - R[i-1] ;
      i = LD$ty;              Opred[i] = S[i] - (R[i-1] + R[i])/2 ;
      i = (LD$ty+1):LD$N ;  Opred[i] = S[i] - R[i] ;
      Opred = .Internal(pmax( na.rm=TRUE, K*q*Opred, LD$eps ) ) ;
      if (LD$Missing$nO > 0) {
        LD$O[LD$Missing$O ] = rlnorm( LD$Missing$nO, parm[LD$i$Omissing], sd=LD$cv )
        LD$O[LD$Missing$O ] = .Internal(pmax(na.rm=TRUE, LD$O[LD$Missing$O ], LD$eps )) ;
        loglik[LD$i$Omissing] = dnorm( parm[LD$i$Omissing], LD$O[LD$Missing$O ], LD$cv*R[LD$Missing$O ] )
      }
      ll_obs = dnorm( log(LD$O), log(Opred), O_sd, TRUE )

      # Removals
    #  for( int i=1; i<nt; i++){
    #    nloglik[2] -= dnorm( R(i-1), Rp(i-1), Rp(i-1)*cv_R, true );
    #  }

    # forecast and monitoring
    #  for( i 1:LD$N ){
    #    ER[i] = R[i] / S[i] ;
      #  B(i) <- S(i)*K
      #  C(i) <- R(i)*K
    #    F[i] = -log(1 - ER[i] ) ; # fishing mortality
    #  }

      LL = sum( ll_obs )  # log likelihood (of the data)
      LP = LL + sum( loglik )  # log iterior
      out = list( LP=LP, Dev=-2*LL, Monitor=c(LP, r, K, q), yhat=S*K, parm=parm )
      return( out  )
    }
    LD$Model <- compiler::cmpfun(LD$Model) #  byte-compiling for more speed .. use RCPP if you want more speed
  }
  return (LD)

}

