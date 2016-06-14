
from numpy import *
from pymc import *
from rpy2.robjects import r
import pdb

# invlogit, Lognormal, deterministic, stochastic, data, Uniform, MCMC, lognormal_like, Gamma, normal_like, poisson_like, potential, uniform_like, observed, Lambda, Normal, AdaptiveMetropolis
# from pymc.Matplot import plot
# from numpy import zeros, ones, array, mean, std, exp, sqrt, pi, log


r('source( "~/.Rprofile" )')
r('p=bio.snowcrab::load.environment()" ) ')
r('res = biomass.summary.db(p=p)')

IOA = array( r('as.matrix(res$B)') )
CAT = array( r('as.matrix(res$L)') )

U = IOA.shape[1]  # number of regions
N = IOA.shape[0]  # no years with data

# b_inits = arrary( )
b0x = array([1., 1., 1./4.])
K0x = array([6., 60., 1.])
r0x = array([1., 1., 1.])
q0x = array([1., 1., 1.])

cv = 0.4 #default 1/cv  -- muliplier for dgamma on inverse scale
er = 0.2  # target exploitation rate
M = 3 # no years for projections
ty = 7  # index of the transition year (2004) between spring and fall surveys
cfa4x = 3 # index of cfa4x
eps = 1e-3  # small non-zero number


# Initialize a few variables/matrices
Kcv=Kmu=Ktau=K=list([None]*U)
rcv=rmu=rtau=r=list([None]*U)
qcv=qmu=qtau=q=list([None]*U)
bcv=b0_mu=btau=list([None]*U)
Otau=Osd=Ocv=list([None]*U)
b=bm=rem=array([[None]*U]*N)
b_inits = array([[0.5]*U]*N)


# Stochastic nodes -- priors + hyperpriors
for j in range(U):
  Kcv[j] = Uniform('Kcv', eps, 1., value=0.1 )
  Kmu[j] = Uniform('Kmu', eps, K0x[j]*1.5, value=0.5 )
  Ktau[j] = Lambda( 'qtau', lambda cv=Kcv[j], mu=Kmu[j] : (cv*mu) **-2 )
  K[j] = Normal('K', Kmu[j], Ktau[j] )

for j in range(U):
  rcv[j] = Uniform('rcv', eps, 1., value=0.1 )
  rmu[j] = Uniform('rmu', eps, 3., value=1 )
  rtau[j] = Lambda( 'rtau', lambda cv=rcv[j], mu=rmu[j] : (cv*mu) **-2 )
  r[j] = Normal('r', rmu[j], rtau[j] )

for j in range(U):
  qcv[j] = Uniform('qcv', eps, 1., value=0.1 )
  qmu[j] = Uniform('qmu', eps, 2., value=1 )
  qtau[j] = Lambda( 'qtau', lambda cv=qcv[j], mu=qmu[j] : (cv*mu) **-2 )
  q[j] = Normal('q', qmu[j], qtau[j] )


# (No) Catch observation model
for j in range(U):
  rem[:,j] = Lambda( 'rem', lambda C=CAT[:,j], KP=K[j] : C/KP )


# Biomass process model
for j in range(U):
  b0_mu[j] =  Uniform('b0_mu', b0x[j]/2., 1.25 )
  bcv[j] = Uniform('bcv', eps, 1., value=0.1 )
  btau = Lambda( 'btau', lambda cv=bcv[j], mu=b0_mu[j] : (cv*mu) **-2 )
  for i in range(1,N):
    if i==0 :
      b[0,j] = Lognormal('b0', mu=b0_mu[j], tau=btau, value=b_inits[0,j])
    else:
      bmean = Lambda("bmean", lambda B=b[i-1,j], R=r[j], C=rem[i-1,j], eps=eps: log(max( B+R*B*(1-B)-C,eps)) )
      b[i,j] = Lognormal('b%i'%i, mu=bmean, tau=btau, value=b_inits[i,j])



# Biomass observation model
#   require a fall / spring surveys correction
#   spring surveys from 1998 to 2003
#   want B to represent the total biomass available in fishing year y
#   when surveys are conducted in fall, Btot(t) = Bsurvey(t) + removals(t)
#   when surveys are conducted in spring, Btot(t) = Bsurvey(t) + removals(t-1)
#   |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
#   corrections for spring surveys in 2003 and earlier for north and south
#   assume BT = total biomass = rem(t-1) + biomass estimate(t; BE)
#   Cfa 4X -- assume total b = rem(t-1) + b estimate(t)
#   ie, similar to a spring fishery but no need for separate q's
#   N- & S-ENS = columns 0 1; CFA 4x column 2

for j in range(U):
  Ocv[j] = Uniform('Ocv', eps, 1., value=0.1 )
  Otau[j] =  Lambda( 'Otau', lambda cv=Ocv[j] : (cv) **-2 ) # assume mu ~ 1
  if j in [0,1]:
    for i in range(N):
      if i == 0:
        pfb = b[i,j] - rem[i,j] # pfb = post fishery biomass
      elif (i >=1 and i < ty):
        pfb = b[i,j] - rem[i-1,j]
      elif i == ty:
        pfb = b[i,j] - (rem[i-1,j]+rem[i,j])/2.
      else:
        pfb = b[i,j] - rem[i,j] 
      bm[i,j] = Lognormal('O%i'%i, log(max(eps, q[j] * K[j] * pfb)), Otau[j], observed=True, value=IOA[i,j] )
  if j in [2]:
    for i in range(N):
      if i == 0:
        pfb = b[i,j] - rem[i,j]
      else:
        pfb = b[i,j] - rem[i-1,j]
      if i <= 6:
        bm[i,j] = Lognormal('O%i'%i, log(max(eps, q[j] * K[j] * pfb)), Otau[j], observed=False )
      else:
        bm[i,j] = Lognormal('O%i'%i, log(max(eps, q[j] * K[j] * pfb)), Otau[j], observed=True, value=IOA[i,j] )



# forecasts



# Parameters on regular scale
K = Lambda('K', lambda k=k: 1/k)  # Lambda converts function to deterministic object

# Variances on regular scale
bsd = Lambda('Sigma2', lambda tau=btau: 1./sqrt(tau) )


# Quantities of interest
MSY=Lambda('MSY', lambda r=r, K=K: r*K/4)

@deterministic
def b1990(b=b[N-1], r=r, C=rem[N-1], k=k): b+r*b*(1-b)-k*C

b1990=Lambda('b1990', lambda b1990=b1990, K=K: b1990*K)


import pymc
import spm1

S = pymc.MCMC(spm1, db='pickle')
S.sample( iter=10000, burn=5000, thin=10 )
pymc.Matplot.plot(S)




# cfa4X have very high errors ... constrain by using SENS estimates
  cv.K[3]  <- cv.K[2]
  cv.p[3]  <- cv.p[2]
  cv.o[3]  <- cv.o[2]


# -------------------
  # monitoring nodes and parameter estimates for output
    for(j in 1:U) {
      Bdrop[j]  <- 1 - step( biomass[N+1,j]-biomass[N,j] ) ; # test if biomass(t) >= biomass(t-1)
      BXMSY[j]  <- 1 - step( biomass[N+1,j]-0.5 ) ; # test if biomass >= BMSY
      BX2MSY[j] <- 1 - step( biomass[N+1,j]-0.25 ) ; # test if biomass >= 1/2 BMSY
      MSY[j]    <- r[j]* K[j] / 4  # maximum height of of the latent productivity (yield)
      BMSY[j]   <- K[j]/2  # biomass at MSY
      FMSY[j]   <- 2 * MSY[j] / K[j] # fishing mortality at MSY
      Fcrash[j] <- 4 * MSY[j] / K[j] # fishing mortality at which the stock will crash
    }


    # -------------------
    # fishing mortality
    # force first year estimate assuming catches in year 0 to be similar to year 1
    for(j in 1:U) {
      for(i in 1:N) {
        F[i,j] <- -log( max(1 - rem[i,j] / biomass[i,j], eps))
      }
      for(i in (N+1):(N+M)) {
        F[i,j] <- -log( max(1 - er * biomass[i-1,j] / biomass[i,j], eps))
      }
    }


    # -------------------
    # annual production
    for(j in 1:U) {
      p[1,j] <- biomass[2,j]- biomass[1,j] + rem[1,j] # approximation
      for (i in 2:(N) ){
        p[i,j] <- (biomass[i+1,j]- biomass[i-1,j])/2 + rem[i,j]  # linear interpolation cancels out the biomass[i,j] term
      }
      for(i in (N+1):(N+M-1)) {
        p[i,j] <- (biomass[i+1,j]- biomass[i-1,j])/2 + er * biomass[i-1,j]   # linear interpolation cancels out the biomass[i,j] term
      }
      p[(N+M),j] <- (biomass[(N+M),j]- biomass[(N+M-1),j]) + er * biomass[(N+M-1),j]   # approximation
    }



    # -------------------
    # recaled estimates

    for(j in 1:U) {
      for(i in 1:(N+M)) {
        B[i,j] <- biomass[i,j]*K[j]
        P[i,j] <- p[i,j]*K[j]
        REM[i,j] <- rem[i,j]*K[j]
        C[i,j] <- catch[i,j]*K[j]
      }
      for(i in 1:M) {
        TAC[i,j] <- catch[N+i,j]*K[j]
      }
    }

}



