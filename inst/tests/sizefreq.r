

  # For details see 
  # http://www.ci.tuwien.ac.at/~gruen/BayesMix/bayesmix-intro.pdf
            
  #          K                          
  # P(yi ) = [SUM]  k Normal(yi | µk , sigma(k)^2 )
  #          k=1
 
  require(coda)
  require(rjags)
  rjags::load.module("dic")
  rjags::load.module("glm")

  loadfunctions( "bayesian" )
  loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
 

  data( fish, package="bayesmix")
  
  fdat = list(
    y = as.vector(fish$fish),
    b0 = 5.625,
    B0inv = 0.1,
    # S0 = 2,
    nu0Half = 10,
    g0Half = 1,
    g0G0Half = 1,
    k = 4,
    N = 256,
    e = c(1, 1, 1, 1)
  )

  
  model.fish = "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(mu[S[i]],tau[S[i]])
    S[i] ~ dcat(eta)
  }
  for (j in 1:k) {
    mu[j] ~ dnorm(b0,B0inv)
    tau[j] ~ dgamma(nu0Half, nu0Half*S0)
  }
  S0 ~ dgamma(g0Half,g0G0Half)
  eta ~ ddirch(e[])
} 
"

# jags.module(names, path)

m = jags.model( file=jags.setup(model.fish), data=fdat, n.chains=3, n.adapt=1000 )  #initiate and adapt steps
# x = jags.samples(m, variable.names=c("mu", "tau", "eta", "S"), n.iter=5000, thin = 10) # sample from posterior

# lapply(x, FUN=function(x) apply(x, MARGIN=c(1), FUN=mean) )
# lapply(x, FUN=function(x) apply(x, MARGIN=c(1), FUN=sd) )

variable.names (m)
coef(m)
dic.samples(m,  n.iter=5000 ) # pDIC

update( m, n.iter=5000)
tomonitor = c("mu", "tau", "eta")

thin =10
niter=5000*thin

yc = coda.samples(m, variable.names=tomonitor, n.iter=niter, thin=thin ) # sample from posterior
# or:
# y = as.mcmc( x )
plot(yc)
densityplot(yc)
autocorr.plot(yc[1])
  autocorr.plot(yc)
    # plot(yc, ask=T)
    # autocorr(yc, lags = c(0, 1, 5, 10, 50), relative=TRUE)

acf( yc[1])
autocorr(yc, lags = c(0, 1, 5, 10, 50), relative=TRUE)

geweke.plot(yc)
gelman.plot(yc[1])


  # sample from the posteriors
   y = jags.samples(m, variable.names=tomonitor, n.iter=niter, thin=thin) # sample from posterior



# ----------------------
# using snowcrab data
 
  det = snowcrab.db( DS ="det.georeferenced" ) 
  
  hvar="cw" #!
  base = 2
  s0 = 10  # ~ ln(10)
  s1 = 175  # ~ ln(160)
  det = det[(det[,hvar]>=s0 & det[,hvar]<=s1),]

  # ii = which(det$sex==0 & det$mat=="1" & as.numeric(substring(as.character(det$trip),6,9))==2008)
  ii = which(det$sex==0 )
  
   
  model.snowcrab = "
model {
  for (i in 1:N) {
    y[i] ~ dnorm( mu[S[i]], tau[S[i]] )
    S[i] ~ dcat( eta )
  }
  for (j in 1:k) {
    mu[j] ~ dnorm( b0, B0inv)
    tau[j] ~ dgamma( nu0Half, nu0Half*S0 )
  }
  S0 ~ dgamma(g0Half,g0G0Half)
  eta ~ ddirch(e[])
} 
"


  fdat = list(
    y = as.vector( det$cw[ii] ),
    b0 = 5.625,
    B0inv = 0.1,
    # S0 = 2,
    nu0Half = 10,
    g0Half = 1,
    g0G0Half = 1,
    k = 4,
    N = 256,
    e = c(1, 1, 1, 1)
  )

  m = jags.model( file=jags.setup( model.snowcrab ), data=fdat, n.chains=3, n.adapt=5000 )  

  # lapply(x, FUN=function(x) apply(x, MARGIN=c(1), FUN=mean) )
  # lapply(x, FUN=function(x) apply(x, MARGIN=c(1), FUN=sd) )

  variable.names (m)
  coef(m)
  dic(m)
  update( m, 1000)

  y = coda.samples(m, variable.names=c("mu", "tau", "eta"), n.iter=5000, thin = 10) # sample from posterior
# or:
# y = as.mcmc( x )
  plot(y)
  densityplot(y)
  autocorr.plot(y[1])

  acf( y[1])
  autocorr(y, lags = c(0, 1, 5, 10, 50), relative=TRUE)

  geweke.plot(y)
  gelman.plot(y[1])




