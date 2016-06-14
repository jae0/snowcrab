
import numpy as np
import pymc as pm
import pdb as pdb

# get data from an R-session
from rpy2.robjects import r
r('source( "~/.Rprofile" )')
r('p=bio.snowcrab::load.environment()" ) ')


r('res = biomass.summary.db(p=p)')
IOA = np.array( r('as.matrix(res$B)') )
CAT = np.array( r('as.matrix(res$L)') )

# additional parameters
U = IOA.shape[1]  # number of regions
N = IOA.shape[0]  # no years with data
b0x = np.array([1., 1., 1./4.])
K0x = np.array([6., 60., 1.])
r0x = np.array([1., 1., 1.])
q0x = np.array([1., 1., 1.])
cv = 0.4 #default 1/cv  -- muliplier for dgamma on inverse scale
er = 0.2  # target exploitation rate
M = 3 # no years for projections
ty = 7  # index of the transition year (2004) between spring and fall surveys
cfa4x = 3 # index of cfa4x
eps = 1e-3  # small non-zero number


# Initialize a few variables/matrices
#Kcv=Kmu=Ktau=K= np.empty(U, dtype=object)
#rcv=rmu=rtau=r= np.empty(U, dtype=object)
#qcv=qmu=qtau=q= np.empty(U, dtype=object)
#bcv=b0_mu=btau= np.empty(U, dtype=object)
#Otau=Osd=Ocv=   np.empty(U, dtype=object)

b_inits = np.array([[0.5]*U]*N)
b = bm = np.array([np.empty(U, dtype=object)]*N)


# Stochastic nodes -- priors + hyperpriors
Kcv = pm.Uniform('Kcv', eps, 1., value=[0.1]*U )
Kmu = pm.Uniform('Kmu', eps, K0x*1.5, value=[0.5]*U )
Ktau = pm.Lambda( 'Ktau', lambda cv=Kcv, mu=Kmu : (cv*mu) **-2 )
K = pm.Normal('K', Kmu, Ktau )

rcv = pm.Uniform('rcv', eps, 1., value=[0.1]*U )
rmu = pm.Uniform('rmu', eps, 3., value=[1]*U )
rtau = pm.Lambda( 'rtau', lambda cv=rcv, mu=rmu : (cv*mu) **-2 )
r = pm.Normal('r', rmu, rtau )

qcv = pm.Uniform('qcv', eps, 1., value=[0.1]*U )
qmu = pm.Uniform('qmu', eps, 2., value=[1]*U )
qtau = pm.Lambda( 'qtau', lambda cv=qcv, mu=qmu : (cv*mu) **-2 )
q = pm.Normal('q', qmu, qtau )

# No catch observation model -- assumed known with no error
rem = pm.Lambda( 'rem', lambda C=CAT, K=K : C/K )

# Biomass process model
b0_mu = pm.Uniform('b0_mu', b0x/2., 1.25, value=b0x )
bcv = pm.Uniform('bcv', eps, 1., value=[0.1]*U )
btau = pm.Lambda( 'btau', lambda cv=bcv, mu=b0_mu : (cv*mu) **-2 )


b[0,:] = pm.Lognormal('b0', mu=b0_mu, tau=btau, value=b_inits[0,:])
for j in range(U):
  for i in range(1,N):
    bmean = pm.Lambda("bmean", lambda B=b[i-1,j], R=r[j], C=rem[i-1,j], eps=eps : 
      np.log( max( B+R*B*(1.-B)-C, eps) ) ) 
    b[i,j] = pm.Lognormal('b%i'%i, mu=bmean, tau=btau[j], value=b_inits[i,j])

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

Ocv = pm.Uniform('Ocv', eps, 1., value=[0.1]*U )
Otau =  pm.Lambda( 'Otau', lambda cv=Ocv : (cv) **-2 ) # assume mu ~ 1

for j in range(U):
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
      bm[i,j] = pm.Lognormal('O%i'%i, np.log(max(eps, q[j] * K[j] * pfb)), Otau[j], observed=True, value=IOA[i,j] )
  if j in [2]:
    for i in range(N):
      if i == 0:
        pfb = b[i,j] - rem[i,j]
      else:
        pfb = b[i,j] - rem[i-1,j]
      if i <= 6:
        bm[i,j] = pm.Lognormal('O%i'%i, np.log(max(eps, q[j] * K[j] * pfb)), Otau[j], observed=False )
      else:
        bm[i,j] = pm.Lognormal('O%i'%i, np.log(max(eps, q[j] * K[j] * pfb)), Otau[j], observed=True, value=IOA[i,j] )

