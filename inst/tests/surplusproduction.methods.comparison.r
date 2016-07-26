
  # surplus production analysis -- comparison of methods

# base data
res = bio.snowcrab::biomass.summary.db()
yrs = as.numeric( rownames( res$B) )
nyr = length(yrs)
O = res$B[ , which(colnames(res$B)=="cfanorth")]  # system/survey index/observations
#O[O>80] = 60
#O[ which( yrs <= 1998 ) ] = 60 # cfa.south.bad.data

R = res$L[, which(colnames(res$L)=="cfanorth")]  # removals/loss



# --------------------------------------
# 1. ML / LaplacesDemon -- shared code

# install_github("lazycrazyowl/LaplacesDemon")
# install_github("lazycrazyowl/LaplacesDemonCpp")
# require(LaplacesDemon)
require(LaplacesDemonCpp)

# data object
LD = list(
  O = O, # observed index of abundance
  removals = R , # removalsches  , assume 20% handling mortality and illegal landings
  er = 0.2,  # target exploitation rate
  N = length( O ) , # no years with data
  M = 5, # no years for projections
  MN = length( O ) + 5,
  ty = which(yrs==2004) ,  # index of the transition year (2004) between spring and fall surveys
  r0= 1,
  K0= 5,
  q0= 1,
  S0= 0.6, # normalised
  cv = 0.5,
  smax =1.25,
  ii = 2:length(O),
  jj = 1:(length(O)-1),
  eps = 1e-6
)

LD$PGF = function(LD) {
  # parameter generating function
  r = rnorm(1, LD$r0, sd=LD$cv )
  K = rlnorm(1, log(LD$K0), sd=LD$cv ) # log scale
  q = rnorm(1, LD$q0, sd=LD$cv )
  q_sd = rhalfcauchy(1, LD$q0 )  
  r_sd = rhalfcauchy(1, LD$r0 ) 
  K_sd = rhalfcauchy(1, LD$K0/5 ) 
  S_sd = rhalfcauchy(1, LD$cv)
  O_sd = rhalfcauchy(1, LD$cv)
  S = rlnorm( LD$N, log(LD$S0), sd=LD$cv )
  S0 = rnorm(1, LD$S0, sd=LD$cv )
  return( c(r, K, q, q_sd, r_sd, K_sd, S_sd, O_sd, S, S0) )
}

# paramater names and initial values
LD$parm.names = as.parm.names( list(
  r = LD$r0,
  K = LD$K0,
  q = LD$q0,
  q_sd = LD$cv*LD$q0,
  r_sd = LD$cv*LD$r0,
  K_sd = LD$cv*LD$K0,
  S_sd = LD$cv,
  O_sd = LD$cv,
  S = rep( LD$S0, LD$N ),
  S0 = LD$S0
))

# index position of paramaters
LD$i = list(
  r=grep("\\<r\\>", LD$parm.names),
  K=grep("\\<K\\>", LD$parm.names),
  q=grep("\\<q\\>", LD$parm.names),
  r_sd=grep("\\<r_sd\\>", LD$parm.names),
  K_sd=grep("\\<K_sd\\>", LD$parm.names),
  q_sd=grep("\\<q_sd\\>", LD$parm.names),
  S_sd=grep("\\<S_sd\\>", LD$parm.names),
  O_sd=grep("\\<O_sd\\>", LD$parm.names),
  S=grep("\\<S\\>", LD$parm.names),
  S0=grep("\\<S0\\>", LD$parm.names)
)

LD$mon.names = c("LP", "r", "K", "q" )
# LD$mon.names = c("LP", "r", "K", "q", paste0("AR",1:(LD$N-1) ) )
# LD$mon.names = c("LP", "r", "K", "q", paste0("S",1:LD$N))


Model = function(parm, LD) {

  Spred = AR = Opred = ER = F = rep(0, LD$N )   # initialize

  # Priors
  loglik = c()
  loglik[LD$parm.names] = 0

  q = parm[LD$i$q] = interval( parm[LD$i$q], LD$eps, Inf )
  r = parm[LD$i$r] = interval( parm[LD$i$r], LD$eps, Inf )
  K = parm[LD$i$K] = interval( parm[LD$i$K], LD$eps, Inf )
  q_sd = parm[LD$i$q_sd] = interval( parm[LD$i$q_sd], LD$eps, Inf )
  r_sd = parm[LD$i$r_sd] = interval( parm[LD$i$r_sd], LD$eps, Inf )
  K_sd = parm[LD$i$K_sd] = interval( parm[LD$i$K_sd], LD$eps, Inf )
  O_sd = parm[LD$i$O_sd] = interval( parm[LD$i$O_sd], LD$eps, Inf )
  S_sd = parm[LD$i$S_sd] = interval( parm[LD$i$S_sd], LD$eps, Inf )
  
  loglik[LD$i$q] = dnorm( q, LD$q0, q_sd,  TRUE ) ;
  loglik[LD$i$r] = dnorm( r, LD$r0, r_sd, TRUE ) ;
  loglik[LD$i$K] = dnorm( K, LD$K0, K_sd, TRUE ) ;
  loglik[LD$i$q_sd] = dhalfcauchy( q_sd, LD$r0/5, log=TRUE)
  loglik[LD$i$r_sd] = dhalfcauchy( r_sd, LD$q0/5, log=TRUE)
  loglik[LD$i$K_sd] = dhalfcauchy( K_sd, LD$K0/5,  log=TRUE)
  loglik[LD$i$S0] = dnorm( parm[LD$i$S0], LD$S0, LD$cv, TRUE ) ;
  loglik[LD$i$O_sd] = dhalfcauchy( O_sd, LD$cv, TRUE );
  loglik[LD$i$S_sd] = dhalfcauchy( S_sd, LD$cv, TRUE );

  S = parm[LD$i$S] = interval( parm[LD$i$S], LD$eps, LD$smax )
  R = LD$removals/K ;
  loglik[LD$i$S] = 0
  for ( i in 1:LD$N) {
    if (i==1) {
      Spred = exp( parm[LD$i$S0] );
    } else {
      AR[i] = 1.0 + r*(1-S[i-1]) ;
      Spred = S[i-1] * AR[i] - R[i-1] ;  # simple logistic
      Spred = .Internal(pmax(na.rm=TRUE, Spred, LD$eps )) ;
    }
    loglik[LD$i$S][i] = dnorm( log(Spred), log(S[i]), S_sd, TRUE );
  }

  # 2. observation model
  Kq = K*q
  i = 1 ;
    Opred[i] = S[i] - R[i] ;
  i = 2:(LD$ty-1)
    Opred[i] = S[i] - R[i-1] ;
  i = LD$ty
    Opred[i] = S[i] - (R[i-1] + R[i])/2 ;
  i = (LD$ty+1):LD$N
    Opred[i] = S[i] - R[i] ;

  Opred = .Internal(pmax( na.rm=TRUE, Kq*Opred, LD$eps ) ) ;
  ll_obs = dnorm( log(O), log(Opred), O_sd, TRUE )

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

  LL = sum( ll_obs )  # no priors ...
  LP = LL + sum( loglik )

  out = list( LP=LP, Dev=-2*LL, Monitor=c(LP, r, K, q), yhat=S*K, parm=parm )

  return( out  )
}

Model <- compiler::cmpfun(Model) #Consider byte-compiling for more speed

ML  = function(...) (Model(...)$Dev / 2)  # i.e. - log likelihood
PML = function(...) (- Model(...)$LP)  #i.e., - log posterior 

# test
Model( parm=LD$PGF(LD), LD=LD)

# maximum likelihood
# Fit.ml = optim( par=LD$PGF(LD), fn=ML, LD=LD, hessian=TRUE, method="CG", control=list(maxit=4000)  )
Fit.ml = optim( par=LD$PGF(LD), fn=ML, LD=LD, hessian=TRUE, control=list(maxit=5000)  )
names(Fit.ml$par ) = LD$parm.names
Fit.ml$par # check convergence
print(sqrt( diag( solve(Fit.ml$hessian) )) ) # asymptotic standard errors


# PML .. better ("penalised ML")
# penalized maximum likelihood .. better but still a little unstable depending on algorithm
#Fit.pml = optim( par=LD$PGF(LD), fn=PML, LD=LD,  method="L-BFGS-B", control=list(maxit=5000), hessian=TRUE )
#Fit.pml = optim( par=LD$PGF(LD), fn=PML, LD=LD,  method="BFGS", control=list(maxit=5000), hessian=TRUE )
Fit.pml = optim( par=LD$PGF(LD), fn=PML, LD=LD,  control=list(maxit=5000), hessian=TRUE )
names(Fit.pml$par ) = LD$parm.names
Fit.pml$par

str(Fit.pml) # check convergence
print(sqrt( diag( solve(Fit.pml$hessian) )) ) # assymptotic standard errors


# Bayesian estimates
Fit = LaplaceApproximation(Model, Data=LD, parm=LD$PGF(LD), Iterations=1000) # quick centering
Fit = LaplaceApproximation(Model, Data=LD, parm=as.initial.values(Fit), Iterations=2000, Method="LBFGS", Stop.Tolerance=1.0E-6)


#Fit = LaplaceApproximation(Model, Data=LD, parm=as.initial.values(Fit), Iterations=5000, Method="TR" )
#Fit = LaplaceApproximation(Model, Data=LD, parm=as.initial.values(Fit), Iterations=5000, Method="SPG" )

# check ACF
Fit = LaplacesDemon(Model, Data=LD, as.initial.values(Fit)=as.initial.values(Fit), Covar=Fit$Covar, Iterations=1000, Status=100, Thinning=1 )
Consort(Fit)
plot(Fit, Data=LD)
# run with correct thinning, etc
Fit = LaplacesDemon(Model, Data=LD, as.initial.values(Fit)=as.initial.values(Fit), Covar=Fit$Covar, Iterations=10000, Status=1000, Thinning=100)
Fit = LaplacesDemon(Model, Data=LD, as.initial.values(Fit)=as.initial.values(Fit), Covar=Fit$Covar, Iterations=30000, Status=1000, Thinning=35, Algorithm="AMWG", Specs=list(B=NULL, n=1000, Periodicity=35 ) )
Fit <- LaplacesDemon(Model, Data=LD, as.initial.values(Fit),
     Covar=Fit$Covar, Iterations=10000, Status=1000, Thinning=105,
     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100,
     n=0, w=1))

# alternate methods
Fit = VariationalBayes(Model, Data=LD, parm=as.initial.values(Fit), Iterations=1000,  Samples=1000, CPUs=5 )
Fit = IterativeQuadrature(Model, Data=LD, parm=as.initial.values(Fit), Iterations=100, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )



# --------------------------------------
# 2. using TMB

tmb.model = paste0( "
#include <TMB.hpp>

template<class Type>
Type posfun(Type x, Type eps, Type &penalty){
    penalty += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
    return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{

// Data
  DATA_VECTOR( O );  // observations
  DATA_VECTOR( removals );  // removals

  // Parameters
  PARAMETER( logr );
  PARAMETER( logK );
  PARAMETER( logq );
  PARAMETER( logS_sd ); // system process error
  PARAMETER( logO_sd ); // system index observation error
//  PARAMETER( logcv_R ); // catch observation error
  PARAMETER_VECTOR( logS );

  Type r=exp(logr);
  Type K=exp(logK);
  Type q=exp(logq);
  Type penalty=0;   // penalty associated with posfun test

  Type S_sd=exp(logS_sd);
  Type O_sd=exp(logO_sd) ;
  vector<Type> S=exp(logS) ;

  // Derived quantities
  int nt = O.size();
  vector<Type> Spred( nt );  //sys state predicted
  vector<Type> Opred( nt );  //observed index predicted
  vector<Type> ER( nt ); // exploitation fraction (harvest rate)
  vector<Type> F (nt) ;
  vector<Type> R (nt) ;

  Type eps=1e-6;
  Type ty=7 ;
  Type Olog=0 ;
  Type Slog=0 ;
  Type S0;


  // Objective function
  vector<Type> nloglik(3);
  nloglik.setZero();

  R = removals/K ;

  // penalties for priors and hyperpriors
  Type qPrior = 1 ;
  Type qSD = 0.2 ;
  nloglik[2] -= dnorm( q, qPrior, qSD, true ) ;

  Type rPrior = 1 ;
  Type rSD = 0.2 ;
  nloglik[2] -= dnorm( r, rPrior, rSD, true ) ;

  Type KPrior = 80 ;
  Type KSD = 10 ;
  nloglik[2] -= dnorm( K, KPrior, KSD, true ) ;

  Type S0Prior = 0.6 ;
  Type S0SD = 0.2 ;
  nloglik[2] -= dnorm( S0, S0Prior, S0SD, true ) ;


  //process model
  for( int i=0; i<nt; i++){
    if (i==0) Spred[0] = S0;
    if (i>0) {
      Spred[i] = S[i-1] * ( 1.0 + r*(1-S[i-1])) - R(i-1) ;  // simple logistic
      // if (model==1) Spred[i] = S[i-1] + r * (1.0 - pow( exp(S[i-1]), theta ) ) - R(i-1) ;  // theta logistic  .. double check parameterizatio ... TODO
    }
    penalty = 0.0 ;
    Spred[i] = posfun( Spred[i], eps, penalty ) ;
    nloglik[0] -= penalty ;
    Slog = log( Spred[i] ) ;
    nloglik[0] -= dnorm( log(S[i]), Slog, S_sd, true ) ;
  }


  //observation model
  for( int i=0; i<nt; i++){
    if ( i==0 )           Opred[i] = K*q*(S[i] - R[i]) ;
    if ( i>0 & i<(ty-1) ) Opred[i] = K*q*(S[i] - R[i-1]) ;
    if ( i==ty )          Opred[i] = K*q*(S[i] - (R[i-1] + R[i])/2 ) ;
    if ( i>ty )           Opred[i] = K*q*(S[i] - R[i]) ;
    penalty = 0.0 ;
    Opred[i] = posfun( Opred[i], eps, penalty ) ;
    nloglik[0] -= penalty ;
    Olog = log( Opred[i] ) ;
    nloglik[1] -= dnorm( log(O[i]), Olog, O_sd, true );
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

  // Reporting
  REPORT( S );
  REPORT( r );
  REPORT( K );
  REPORT( q );
  REPORT( ER );
  REPORT( F );
  REPORT( R );
  REPORT( Opred );
  REPORT( O_sd );
  REPORT( S_sd );

  ADREPORT( S );
  ADREPORT(nloglik);

  return nllik;
}
")



Data = list( O=O, removals=R )
N =  length(O)
Params = list( logr=log(1), logK=log(5),  logq=log(1),
               logS_sd=log(0.5), logO_sd=log(0.5), logS=log( rep(0.6, N)) )


setwd("/tmp")
wd = getwd()
fn = tempfile(pattern="tmb_", tmpdir=wd, fileext=".cpp" )
print(fn)
cat( tmb.model, file=fn)
TMB::compile( fn )
shlib_fn = gsub( ".cpp$", "", basename(fn) )
# dyn.unload(TMB::dynlib( shlib_fn ) )
dyn.load( TMB::dynlib( shlib_fn ) )
Obj = TMB::MakeADFun( data=Data, parameters=Params, random=c("logS" ) )
Obj$fn( Obj$par ) # make sure program will run properlyL-BFGS-B
Obj$gr( Obj$par ) # make sure program will run properlyL-BFGS-B

Opt = optim( par=Obj$par, fn=Obj$fn, gr=Obj$gr,
  method ="Nelder-Mead",
  control=list(trace=1, eval.max=1e5, iter.max=1e5))
plot( Obj$report()$S )

SD = TMB::sdreport( Obj, bias.correct=TRUE)
(SD)


        Estimate Std. Error
logr     0.006065007  0.1967278  -- 1.0
logK     4.301091004  0.1355165  -- 73.8
logq    -0.110785527  0.2213573  -- 0.895


  # --------------------------------------
  # 3. Bayesian via jags/mcmc


logistic.model.simple.jags = paste0("
model {
  cv.K ~ dunif( eps, 1 )
  cv.q ~ dunif( eps, 1 )
  cv.r ~ dunif( eps, 1 )
  cv.p ~ dunif( eps, 1 )
  cv.o ~ dunif( eps, 1 )

  # priors of key stochastic nodes for estimation
  q ~ dnorm( q0, pow( q0 * cv.q , -2 ) ) T(eps,)
  r ~ dnorm( r0, pow( r0 * cv.r , -2 ) ) T(eps,)
  K ~ dnorm( K0, pow( K0 * cv.K , -2 ) ) T(eps,)

  # removals... no observation model, standardized to K
  for (i in 1:N){
    R[i] <- removals[i]/K
  }

  # biomass observation model
  # spring surveys from 1998 to 2003

  Kq <- K * q # multiply here to have fewer operations below

  O[1] ~ dlnorm( log( max( Kq*(S[1] - R[1]) , eps)), pow( cv.o , -2 ) )  # approximation
  for (i in 2:(ty-1)) {
    O[i] ~ dlnorm( log( max( Kq*(S[i]- R[(i-1)]), eps)), pow( cv.o , -2 ) )  ;
  }
  # transition year
  O[ty] ~ dlnorm( log( max( Kq*(S[ty] - (R[(ty-1)] + R[ty] )/2 ), eps)), pow( cv.o , -2 ) ) ;  # approximation
  # fall surveys
  for (i in (ty+1):N) {
    O[i] ~ dlnorm( log( max( Kq*(S[i] - R[i]), eps)), pow( cv.o , -2 ) ) ;
  }

  # S - system state: process model
  # b0 is the starting b which is infered as well
  S[1] ~ dlnorm( log(S0), pow( S0*cv.p, -2) ) T(eps, smax ) ;  # S at first year

  # using observed abundance index CV's -- these are already on the log scale
  for(i in 2:(N+M)) {
    S[i] ~ dlnorm( log( max(S[i-1]*( 1 + r*(1-S[i-1])) - R[i-1] , eps)), pow( cv.p , -2 ) ) T(eps, smax) ;
  }

  # forecasted / predicted TACs at target er
  for(i in 1:M) {
    R[N+i] <- er*S[N+i-1]
  }

# recaled estimates
  for(i in 1:(N+M)) {
    B[i] <- S[i]*K
    C[i] <- R[i]*K
    F[i] <- -log( max(1 - R[i] / S[i], eps)) # fishing mortality
  }

  # monitoring nodes and parameter estimates for output
#  MSY    <- r* K / 4  # maximum height of of the latent productivity (yield)
#  BMSY   <- K/2  # S at MSY

#  FMSY   <- 2 * MSY / K # fishing mortality at MSY

}
")

fn = tempfile(pattern="jags_", tmpdir=getwd(), fileext=".jags" )
cat( logistic.model.simple.jags, file=fn)


Data = list(
  O = O, # observed index of abundance
  removals = R , # removalsches  , assume 20% handling mortality and illegal landings
  er = 0.2,  # target exploitation rate
  N = length( O ) , # no years with data
  M = 5, # no years for projections
  MN = length( O ) + 5,
  ty = which(yrs==2004) ,  # index of the transition year (2004) between spring and fall surveys
  r0=1,
  K0=log( max(O, na.rm=TRUE)),
  q0=1,
  S0=0.5, # normalised
  smax =1.25,
  Kmin = median(O, na.rm=TRUE),
  Kmax = max(O, na.rm=TRUE) * 2,
  eps = 1e-6  # small non-zero number
)


m = rjags::jags.model( file=fn, data=Data, n.chains=3, n.adapt=1000 )

tomonitor =  c( "r", "K", "q", "S", "MSY", "BMSY", "FMSY", "F", "C", "B" )
tomonitor = intersect( variable.names (m), tomonitor )

# MCMC/Gibbs parameters
n.adapt = 1000 # burn-in  .. 4000 is enough for the full model but in case ...
n.iter = 1000
n.chains = 3
n.thin = 200 # use of uniform distributions causes high autocorrelations ?
n.iter.final = n.iter * n.thin


y = rjags::jags.samples(m, variable.names=tomonitor, n.iter=n.iter.final, thin=n.thin) # sample from posterior

apply(y$r, 1, mean, na.rm = T)
apply(y$q, 1, mean, na.rm = T)
apply(y$K, 1, mean, na.rm = T)

> apply(y$r, 1, mean, na.rm = T)
[1] 0.7689612
> apply(y$q, 1, mean, na.rm = T)
[1] 1.015021
> apply(y$K, 1, mean, na.rm = T)
[1] 78.07672


hist(y$r, "fd")
hist(y$q, "fd")
hist(y$K, "fd")

yrsall = c(yrs, max(yrs)+(1:Data$M))
plot(apply(y$B, 1, mean, na.rm = T)~yrsall, type="b", ylim=quantile( c( y$B, Data$O ), probs=c(0.05, 0.95),na.rm=TRUE) )
plot(apply(y$F, 1, mean, na.rm = T)~yrsall, type="b" )

points( O~yrs, Data, col="red", pch=20)

# ----------------
# convergence testing -- by 1000 to 1500 convergence observed by Gelman shrink factor diagnostic

ydic = rjags::dic.samples(m, n.iter=n.iter ) # pDIC

y = rjags::jags.samples(m, variable.names=tomonitor, n.iter=60000, thin=1000 )

coda::gelman.plot(y[["r"]])
coda::gelman.plot(y[["K"]])
coda::gelman.plot(y[["q"]])  # about 6-8000 runs required to converge


# ------------------
# determine autocorrelation thinning
yc = rjags::coda.samples(m, variable.names=c("K", "r", "q"), n.iter=20000, thin=1) # sample from posterior
coda::autocorr.plot(yc)   # about 10 to 20 required
# plot(yc, ask=T)

# coda::geweke.plot(yc[["r"]])

# coda::autocorr(yc, lags = c(0, 1, 5, 10, 50), relative=TRUE)



# --------------------------------------
# 4. Bayesian via Stan --- not working negative values are slipping in ..
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_code = paste0('

data {
  int N ; // number of years with data
  int M ; // number of years to project
  // int MN ; // M+N
  int ty ; // year index of transition
  real r0  ; //  prior mean
  real K0  ; //  prior mean
  real q0  ; //  prior mean
  real S0 ; //  prior mean
  vector[N] O; //  observed index of abundance
  vector[N] removals; //  removal/catch
  real er ; // target exploitation rate
  real eps ; // level of numeric tolerance
  real smax; // max value of S
  real Kmin;
  real Kmax;
  real cv;
}

transformed data {
  vector[N] logO ;
  logO = log(O);
}

parameters {
  real K;
  real r;
  real q;
  real K0sd;
  real r0sd;
  real q0sd;
  vector[N] S;
  vector[N] Ssd;
  vector[N] Osd;
}

transformed parameters {
}

model {
  vector[N] Spred ;
  vector[N] Opred;

  q0sd ~ uniform( eps, cv ) ;
  r0sd ~ uniform( eps, cv ) ;
  K0sd ~ uniform( eps, cv*K0 ) ;

  Ssd  ~ uniform( eps, cv )  ;
  Osd  ~ uniform( eps, cv*K0 )  ;

  // # priors of key stochastic nodes for estimation
  q ~ normal( exp(q0), q0sd ) T[eps,];
  r ~ normal( exp(r0), r0sd ) T[eps,];
  K ~ normal( exp(K0), K0sd ) T[Kmin,Kmax];

  for (i in 1:N) {
    if (i==1)     Spred[i] ~ normal( exp(S0), cv ) T[eps,smax]; // guess, S0 is log scale .. exp ensures positive values
    if (i>1)      Spred[i] = S[i-1]*( 1 + r*(1-S[i-1])) - removals[i-1]/K ;
    Spred[i] = fmax( eps, Spred[i] );
    S[i] ~ lognormal( Spred[i], Ssd[i] ) T[eps,smax];
  }

  // observation model
  for (i in 1:N) {
    if (i==1)               Opred[i] = K * q * (S[i] - removals[i]/K)   ; // # approximation
    if ((i>=2) && (i <= (ty-1))) Opred[i] = K * q * (S[i] - removals[i-1]/K ) ; // spring surveys from 1998 to 2003
    if (i==ty)              Opred[i] = K * q * (S[i] - (removals[i-1] + removals[i] )/(2*K) );
    if (i>ty)               Opred[i] = K * q * (S[i] - removals[i]/K) ;
    Opred[i] = fmax( eps, Opred[i]  ); // expected biomass after fishing
    O[i] ~ normal( Opred[i], Osd[i] ) T[eps,Kmax];
  }

  // forecast
  //for (i in (N+1):M) {
  //  Spred = S[i-1]*( 1 + r*(1-S[i-1])) - S[i-1]*er ; // forecast
  //  print("i:",i  ) ;
    // Spred = fmax( eps, Spred );
  //  print("Spred:", Spred ) ;
  //  Ssd[i] ~ cauchy( 0, 0.25 );  // half-cauchy
  //  S[i] ~ normal( ( Spred), Ssd[i]  ) T[eps,smax];
   // print( "S=", S[i] );
  // }
}

/*
generated quantities {
  vector[MN] REM; //  removal/catch again as previous in model namespace only
  vector[MN] B;
  vector[MN] C;
  vector[MN] F;
  vector[MN] ER;
  real MSY;
  real FMSY;
  real BMSY;

  // rescale to original scale
  for (i in 1:MN){
    if (i <= N) REM[i] = removals[i]/K;  // redo as the scaling was done locally above
    if (i > N) REM[i] = er * S[i-1] ;  // forecast
    B[i] = K*S[i];
    C[i] = REM[i]*K ;
    ER[i] = REM[i] / S[i] ;
    if ((1.0 - ER[i]) > eps) {
      F[i]  = -log(eps) ; //# fishing mortality
    } else {
      F[i]  = - log( 1.0 - ER[i] ) ; //# fishing mortality
    }
  }

  // # monitoring nodes and parameter estimates for output
  MSY  = r* K / 4.0  ; //  maximum height of of the latent productivity (yield)
  BMSY = K/2.0  ; // S at MSY
  FMSY = 2.0 * MSY / K ; // fishing mortality at MSY

}
*/

')


# data object
Data = list(
  O = O, # observed index of abundance
  removals = R , # removalsches  , assume 20% handling mortality and illegal landings
  er = 0.2,  # target exploitation rate
  N = length( O ) , # no years with data
  M = 5, # no years for projections
  MN = length( O ) + 5,
  ty = which(yrs==2004) ,  # index of the transition year (2004) between spring and fall surveys
  r0= log(1),
  K0= log( max(O, na.rm=TRUE)),
  q0= log(1),
  S0= log(0.6), # normalised
  cv = 0.5,
  smax =1.25,
  Kmin = median(O, na.rm=TRUE),
  Kmax = max(O, na.rm=TRUE) * 2,
  ii = 2:length(O),
  jj = 1:(length(O)-1),
  eps = 1e-6mO
)


initfn = function() {
  list(
    K=Data$K0, r=Data$r0, q=Data$q0, K0sd=Data$cv, r0sd=Data$cv, q0sd=Data$cv,
    S=rep(Data$S0,Data$N), Ssd=rep(Data$cv,Data$N), Osd=rep(Data$cv,Data$N)
  )
}
# test --- some exception error ... 
fit <- stan(model_code=stan_code, data=Data, iter=1000, chains=1,  seed=1, verbose=TRUE, init=initfn)

fit <- stan(model_code=stan_code, data=Data, iter=1000, chains=3, init=initfn, verbose=FALSE)

print(fit)
plot(fit)
pairs(fit, pars = c("K", "r", "q"))
pairs(fit, pars = c("S", "Ssd"))

la <- extract(fit, permuted = TRUE) # return a list of arrays
mu <- la$mu

### return an array of three dimensions: iterations, chains, parameters
a <- extract(fit, permuted = FALSE)

### use S3 functions as.array (or as.matrix) on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)







# --------------------------------------
# 6. ML via Julia Optim

--to do



# --------------------------------------
# 7. Bayesian via Julia / Mamba

--to do









