

require(LaplacesDemonCpp)

  
if (!exists("current.year")) current.year=year(Sys.Date())
p = bio.snowcrab::load.environment( year.assessment=current.year)

debug.region="cfa4x"
debug.region="cfasouth"
debug.region="cfanorth" 

sb = surplusproduction.db( DS="LaplacesDemon.debug", sourcedata="nosa", debug.region=debug.region) 
sb = surplus.production.simple.laplacesdemon.setup( sb )   # set up the model



# 2. maximum likelihood solution
f.ml = optim( par=sb$PGF(sb), fn=sb$Model.ML, Data=sb, method="BFGS", control=list(maxit=5000, trace=0), hessian=TRUE  )
names(f.ml$par ) = sb$parm.names
#print(sqrt( diag( solve(f.ml$hessian) )) ) # assymptotic standard errors
(f.ml$par)

# 3. penalized maximum likelihood .. better but still a little unstable depending on algorithm
f.pml = optim( par=sb$PGF(sb), fn=sb$Model.PML, Data=sb, method="BFGS", control=list(maxit=5000, trace=0), hessian=TRUE )
names(f.pml$par ) = sb$parm.names
(f.pml$par)

#print(sqrt( diag( solve(f.pml$hessian) )) ) # assymptotic standard errors



f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=1000, Status=100, Thinning=10) # 

f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=10) # 

f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=as.initial.values(f), Iterations=5000, Status=10, Thinning=20, Method="NUTS") #


Initial.Values <- as.initial.values(f)
f <- LaplacesDemon(sb$Model, Data=sb, Initial.Values,
     Covar=f$Covar, Iterations=20000, Status=1000, Thinning=1000,
     Algorithm="CHARM", Specs=NULL)


f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=as.initial.values(f), Iterations=5000, Status=100, Thinning=25, Covar=f$Covar)  



f = LaplaceApproximation(sb$Model, Data=sb, parm=sb$PGF(sb), Method="Roptim", method="BFGS", Stop.Tolerance=1e-9, Iterations = 10000  )
f

Consort(f)
plot(f, Data=sb)

PosteriorChecks(f)


f0 = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=1000, Status=100, Thinning=1)


f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="HAR", Iterations=1000  ) 

f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="SR1", Iterations=1000  ) 

f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="TR", Iterations=10000  ) 


f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="BFGS", Iterations=1000  ) 


f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="PSO", Iterations=10000  ) 


f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="SPG", Iterations=10000  ) 




# 4. quick solution: acts as "burn-in" .. do a few times in case solution is unstable
f = LaplaceApproximation(sb$Model, Data=sb, parm=sb$PGF(sb), Iterations=100 ) 
f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Iterations=1000 ) 
f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="BFGS", Stop.Tolerance=1e-8, Iterations=1000  ) 


f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=5000, Status=10, Thinning=10, Algorithm="NUTS", Covar=f$Covar, Specs=list(A=100, delta=0.6, epsilon=NULL, Lmax=Inf))  # A=burnin, delta=target acceptance rate

f = LaplacesDemon.hpc(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=10, Status=1, Thinning=1, Algorithm="NUTS", Covar=f$Covar, Specs=list(A=100, delta=0.6, epsilon=NULL, Lmax=Inf))  # A=burnin, delta=target acceptance rate

f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=1000, Status=100, Thinning=1)

f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="Roptim", Stop.Tolerance=1e-8, Iterations = 5000, method="Nelder-Mead", control=list(maxit=5, reltol=1e-9 ) )
f

f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="Roptim", Stop.Tolerance=1e-8, Iterations = 1000, method="L-BFGS-B", control=list(maxit=10, reltol=1e-10 ) )



f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f)  ) # fast spin up of paramters
f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="BFGS", Stop.Tolerance=1e-8  ) 


# NOTES: NM is slow
# algorithms = c( "TR", "LM", "CG", "NR", "NM", "SPG", "LBFGS", "BFGS", "PSO", "SR1", "HAR", "DFP", "BHHH", "HJ", "Rprop" ) # available
# algorithms = c( "CG", "NM", "SPG", "LBFGS", "BFGS", "PSO", "SR1", "HAR", "DFP", "Rprop" ) # reliably working
f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=500, Status=100, Thinning=1)

algorithms = c( "CG", "LBFGS", "BFGS", "PSO", "SR1", "HAR", "DFP", "Rprop" ) # reliably working
for (a in algorithms) {
  print(a)
  ft = try( LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method=a, Stop.Tolerance=1e-9, Iterations=5000 ) )
  if (! class(ft) %in% "try-error" ) if (ft$Converged) if( ft$LP.Final > f$LP.Final )  {
    f = ft
    f$Method = a 
  }
}
str(f)
plot(f, Data=sb)

# MCMC ... look at ACF
f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=as.initial.values(f), Covar=f$Covar, 
  Iterations=1000, Status=100, Thinning=1, Algorithm="CHARM" )
Consort(f)
plot(f, Data=sb)

 PosteriorChecks(f)
     #caterpillar.plot(f, Parms="beta")
     #plot(f, Data, PDF=FALSE)
     #Pred <- predict(f, Model, Data, CPUs=1)

  #summary(Pred, Discrep="Chi-Square")
     #plot(Pred, Style="Covariates", Data=Data)
     #plot(Pred, Style="Density", Rows=1:9)
     #plot(Pred, Style="Fitted")
     #plot(Pred, Style="Jarque-Bera")
     #plot(Pred, Style="Predictive Quantiles")
     #plot(Pred, Style="Residual Density")
     #plot(Pred, Style="Residuals")
     #Levene.Test(Pred)
     #Importance(f, Model, Data, Discrep="Chi-Square")

# MCMC: run with appropriate thinning and other options:
f = LaplacesDemon(sb$Model, Data=sb, as.initial.values(f),
  Covar=NULL, Iterations=10000, Status=1000, Thinning=1000, Algorithm="CHARM", Specs=NULL)

f = LaplacesDemon(sb$Model, Data=sb, as.initial.values(f),
  Covar=f$Covar , Iterations=50000, Status=10204, Thinning=1000, Algorithm="CHARM", Specs=list(alpha.star=0.44))


f = LaplacesDemon(sb$Model, Data=sb, 
  Initial.Values=as.initial.values(f), Covar=f$Covar, Iterations=10000, Status=1000, Thinning=100)

f = LaplacesDemon(sb$Model, Data=sb, 
  Initial.Values=as.initial.values(f), Covar=f$Covar, Iterations=30000, Status=1000, Thinning=35, 
  Algorithm="AMWG", Specs=list(B=NULL, n=1000, Periodicity=35 ) )

f <- LaplacesDemon(sb$Model, Data=sb,
  Initial.Values=as.initial.values(f), Covar=f$Covar, Iterations=10000, Status=1000, Thinning=105, 
  Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))

# medium speed .. increase Iterations til convergence
f = VariationalBayes(sb$Model, Data=sb,  parm=as.initial.values(f), Iterations=100,  Samples=10, CPUs=5 )
f = VariationalBayes(sb$Model, Data=sb,  parm=as.initial.values(f), Iterations=5000,  Samples=100, CPUs=5 ,Stop.Tolerance=1e-9)

# slow
f = IterativeQuadrature(sb$Model, Data=sb, parm=as.initial.values(f), Iterations=10, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )





--------------

sb = surplusproduction.db( DS="jags.2015", sourcedata="nosa" ) 
# sb = surplusproduction.db( DS="jags.2015", sourcedata="default" ) 
