
using Distributions
using Optim


# data
O = [ 60.00000 60.00000 60.00000 77.53333 62.65507 31.18533 30.62145 41.31591 32.66703 40.07025 39.74323 51.51829 33.07368 32.29060 35.07408 42.91545 ]

removals = [ 1.558000  2.700000  8.701000  9.048000  8.891000  8.836000  8.095283  6.412215  4.489686 4.946139  8.265875 10.760470 12.760896 12.160676 11.717945 11.351338 ]


# constants
r0 =  1
K0 = 75
q0 = 1
S0 = 0.6
er = 0.2
N = length(O)
M =  5
MN = M + N
ty = 7
eps = 1e-06

qSD = 0.2 ;
rSD = 0.2 ;
KSD = 10 ;


Spred = rand(Uniform(), N)
Osd   = rand(Uniform(), N)
Ssd   = rand(Uniform(), N)


# params
pm = [ K0; r0; q0; S0; Spred; Osd; Ssd ]

# indices of params
iK = 1
ir = 2
iq = 3
iS0 = 4
iS = 5 + [1:N]
iOsd = maximum(iS) + [1:N]
iSsd = maximum(iOsd) + [1:N]


function bd(pm::Vector)

  # unpack and transform params
  r=exp(pm[ir]);
  K=exp(pm[iK]);
  q=exp(pm[iq]);
  S_sd=exp(pm[iSsd]);
  O_sd=exp(pm[iOsd] ;
  S=exp(pm[iS]) ;

  R = removals/K ;
  // penalties for priors and hyperpriors

  // priors of key stochastic nodes for estimation
  nloglik[2] -= dnorm( q, qPrior, qSD, true ) ;
  nloglik[2] -= dnorm( r, rPrior, rSD, true ) ;
  nloglik[2] -= dnorm( K, KPrior, KSD, true ) ;
  nloglik[2] -= dnorm( S0, S0Prior, S0SD, true ) ;

  //process model
  for( int i=0; i<nt; i++){
    if (i==0) Spred[0] = S0;
    if (i>0) {
      Spred[i] = S[i-1] * ( 1.0 + r*(1-S[i-1])) - R[i-1] ;  // simple logistic
      // Spred[i] = S[i-1] * (1.0 + r*( 1 - pow( exp(S[i-1]), theta ) )) - R[i-1] ;  // theta logistic  .. double check parameterizatio ... TODO
    }
    Spred[i] = posfun( Spred[i], eps, penalty ) ;
    nloglik[0] -= penalty ;
    Slog = log( Spred[i] ) ;
    if( !isNA(Slog)) nloglik[0] -= dnorm( log(S[i]), Slog, S_sd[i], true ) ;
  }

  //observation model
  for( int i=0; i<nt; i++){
    if ( i==0 )           Opred[i] = K*q*(S[i] - R[i]) ;
    if ( i>0 & i<(ty-1) ) Opred[i] = K*q*(S[i] - R[i-1]) ;
    if ( i==ty )          Opred[i] = K*q*(S[i] - (R[i-1] + R[i])/2 ) ;
    if ( i>ty )           Opred[i] = K*q*(S[i] - R[i]) ;
    Opred[i] = posfun( Opred[i], eps, penalty ) ;
    nloglik[0] -= penalty ;
    Olog = log( Opred[i] ) ;
    if( !isNA(Olog)) nloglik[1] -= dnorm( log(O[i]), Olog, O_sd[i], true );
  }



  for( int i=0; i<nt; i++){
    ER[i] = R[i] / S[i] ;
  //  B(i) <- S(i)*K
  //  C(i) <- R(i)*K
    F[i] = -log(1 - ER[i] ) ; // fishing mortality
  }

  // forecast

  // Removals
//  for( int i=1; i<nt; i++){
//    nloglik[2] -= dnorm( R(i-1), Rp(i-1), Rp(i-1)*cv_R, true );
//  }

  // Total likelihood
  Type nllik = nloglik.sum();

}



end


res = optimize( bd, pm, autodiff=true )

