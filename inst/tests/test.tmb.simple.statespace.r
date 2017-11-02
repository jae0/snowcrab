
# the following is based upon work posted by
# James Thorsen for his course:
# FSH 507 "Spatio-temporal models for ecologists" in
# Spring quarter 2016 at University of Washington
# https://github.com/James-Thorson/2016_Spatio-temporal_models

if (0) {
  devtools::install_github("kaskr/adcomp/TMB")
  install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/testing")
  install.packages( "RandomFields" )
  install.packages( "RANN" )
}


# sim dimensions
nx = 200
ny = 100
nt = 10
nsample = 200
si = sample( 1:(nx*ny), nsample )


locs <- list( x= seq(1, nx), y=seq(1, ny) )
locsGrid <- expand.grid( x=locs$x, y=locs$y )
# unconditional simulation
require(gstat)
modG <- gstat::gstat(formula = z~1, locations = ~x+y, data=locsGrid, dummy=TRUE, beta = 0,
    model = gstat::vgm(psill=0.8, nugget=0.2, model="Mat", range=nx/10, kappa=0.5), nmax = 20)
rfG <- predict(modG, newdata=locsGrid, nsim = 1)

z = rfG[si,3]
locs_xy = rfG[ si, c("x","y") ]

sp::gridded(rfG) = ~x+y
sp::spplot(rfG)

reqiure(geoR)
rfGeoR = geoR::grf( nsim=1, cov.model="matern", cov.pars=c(sigmasq=0.8, phi=0.75), kappa=1)

RandomFields methods ...

require(lbm)
# est cov params
vg = lbm::lbm_variogram( locs_xy, z, methods="gstat" )
vg = lbm::lbm_variogram( locs_xy, z, methods="geoR" )

# prep data for spTimer





tmb.model = paste0( "

#include <TMB.hpp>

// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{

// Data
  DATA_VECTOR( y_t );
  DATA_VECTOR( c_t );
  DATA_VECTOR( ysd_t );
  DATA_VECTOR( penalties_z );

  // Parameters
  PARAMETER( log_r );
  PARAMETER( log_k );
  PARAMETER( log_q );
  PARAMETER( log_sigmap );
  PARAMETER( log_sigmam );
  PARAMETER( log_sigmac );
  PARAMETER_VECTOR( log_x_t );
  PARAMETER_VECTOR( logit_exploit_t );

  // Derived quantities
  int n_years = y_t.size();
  Type r = exp(log_r);
  Type k = exp(log_k);
  Type q = exp(log_q);
  vector<Type> xpred_t( n_years );
  vector<Type> cpred_t( n_years );
  vector<Type> sigmam_t( n_years );
  vector<Type> exploit_t = invlogit( logit_exploit_t );
  vector<Type> x_t = exp( log_x_t );
  vector<Type> Depletion_t = x_t / k;
  cpred_t.setZero();

  // Objective funcction
  vector<Type> jnll_comp(5);
  jnll_comp.setZero();

  // Reconstruct time series
  //cpred_t(0) = x_t(0);
  for( int t=1; t<n_years; t++){
    cpred_t(t-1) = x_t(t-1) * exploit_t(t-1);
    xpred_t(t) = (x_t(t-1)-cpred_t(t-1)) * (1.0 + r*(1-(x_t(t-1)-cpred_t(t-1))/k));
    jnll_comp(0) -= dnorm( x_t(t), xpred_t(t), xpred_t(t)*exp(log_sigmap), true );
  }

  // Probability of data
  for( int t=0; t<n_years; t++){
    sigmam_t(t) = pow(( exp(2.0*log_sigmam) + pow(ysd_t(t),2.0) ), 0.5);
    jnll_comp(1) -= dnorm( y_t(t), q*x_t(t), q*x_t(t)*sigmam_t(t), true );
  }

  // Penalty on catches
  for( int t=1; t<n_years; t++){
    jnll_comp(2) -= dnorm( c_t(t-1), cpred_t(t-1), cpred_t(t-1)*exp(log_sigmac), true );
  }

  // Penalty on starting biomass
  jnll_comp(3) -= dnorm( log(Depletion_t(0)), Type(0.0), exp(penalties_z(0)), true );

  // Penalty on interannual variability in exploitation ratio
  for( int t=1; t<(n_years-1); t++){
    jnll_comp(3) -= dnorm( logit_exploit_t(t), logit_exploit_t(t-1), exp(penalties_z(1)), true );
  }

  // Jacobian for log-random effects
  //jnll_comp(4) -= log_x_t.sum();

  // Total likelihood
  Type jnll = jnll_comp.sum();

  // Reporting
  REPORT( log_x_t );
  REPORT( n_years );
  REPORT( cpred_t );
  REPORT( x_t );
  REPORT( xpred_t );
  REPORT( jnll_comp );
  REPORT( r );
  REPORT( k );
  REPORT( q );
  REPORT( exploit_t );
  REPORT( Depletion_t );
  REPORT( log_sigmam );
  REPORT( log_sigmap );

  ADREPORT( log_x_t );
  ADREPORT( x_t );
  ADREPORT( Depletion_t );

  return jnll;
}
")


fn = tempfile(pattern="tmb_", tmpdir=getwd(), fileext=".cpp" )
cat( tmb.model, file=fn)




library( TMB )
TMB::compile( fn )
dyn.load( dynlib( strsplit(basename(fn), "\\.")[[1]][1] ) )
Obj = MakeADFun( data=Data, parameters=Params, random=Random )
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list(trace=1, eval.max=1e4, iter.max=1e4))
Opt[["final_diagnostics"]] = data.frame( "Name"=names(Obj$par), "final_gradient"=Obj$gr(Opt$par))
Report = Obj$report()
unlist( Report[c('Range','SigmaO','SigmaE')] )
SD = sdreport( Obj, bias.correct=TRUE)




setwd( "C:/Users/James.Thorson/Desktop/Project_git/state_space_production_model/" )
library( TMB )

# Parameters
r = 0.5
k = 1000
sigmap = 0.05
sigmam = 0.1
n_years = 100
n_rep = 100
overdispersion = 2  # 1=No overdispersion

#
Date = Sys.Date()
  DateDir = paste0(getwd(),"/",Date,"/")
  dir.create(DateDir)

# Compile results
Results = array(NA, dim=c(2,n_rep,3), dimnames=list(c("without","with"),paste0("Rep_",1:n_rep),c("k","RE_final_depletion","sigmam")))

# Loop
for(repI in 1:n_rep){
  set.seed(repI)

  # exploitation fraction
  exploit_max = r * 0.6
  exploit_t = c( rep(0.01,ceiling(n_years/4)),seq(0.01,exploit_max,length=ceiling(n_years/4)),rep(exploit_max,floor(n_years/4)),seq(exploit_max,r/4,length=floor(n_years/4)) )

  # Simulate
  catch_t = n_t = rep(NA,n_years-1)
  n_t[1] = rnorm(1, k, sigmap)
  for(t in 2:n_years){
    catch_t[t-1] = n_t[t-1] * exploit_t[t-1]
    n_t[t] = (n_t[t-1]-catch_t[t-1]) * (1 + r*(1-(n_t[t-1]-catch_t[t-1])/k)) + rnorm(1, mean=0, sd=sigmap*n_t[t-1])
    if( n_t[t]< k/100 ) n_t[t] = k/100
  }
  nobs_t = n_t + rnorm(n_years, mean=0, sd=sigmam*n_t)

  # Visualize
  png( paste0(DateDir,"Rep_",repI,"_true_dynamics.png"), width=10, height=7.5, res=200, units="in")
    par(mfrow=c(1,2))
    matplot( cbind(n_t,nobs_t), type="l", ylim=c(0,k*1.5))
    plot(catch_t)
  dev.off()

  #################
  # TMB prep
  #################

  # Compile
  Version_statespace = "statespace_v1"
  compile( paste0(Version_statespace,".cpp") )

  # Fit the model
  Data = list("y_t"=nobs_t, "c_t"=catch_t, "ysd_t"=rep(sigmam/overdispersion,n_years), "penalties_z"=c(log(0.1),log(1)) )
  Params = list("log_r"=log(r), "log_k"=log(max(nobs_t)), "log_q"=log(1), "log_sigmap"=log(1), "log_sigmam"=log(1e-10), "log_sigmac"=log(0.01), "log_x_t"=log(nobs_t), "logit_exploit_t"=rep(0,n_years-1) )
  Random = c("log_x_t", "logit_exploit_t")

  # Save plot
  Plot_Fn = function( report, sdsummary, dir, name, ... ){
    png( paste0(dir,"/",name), width=10, height=7.5, res=200, units="in")
      par( mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
      # Biomass
      matplot( cbind(n_t,report$x_t), type="l", ylim=c(0,k*1.5), ylab="", xlab="year", main="Biomass")
      Mat = cbind( report$x_t, sdsummary[which(rownames(sdsummary)=="x_t"),"Std. Error"])
      polygon( x=c(1:n_years,n_years:1), y=c(Mat[,1]+Mat[,2],rev(Mat[,1]-Mat[,2])), col=rgb(1,0,0,0.2), border=NA)
      # Depletion
      matplot( cbind(n_t/k,report$Depletion_t), type="l", ylim=c(0,1.5), ylab="", xlab="year", main="Relative biomass")
      Mat = cbind( report$Depletion_t, sdsummary[which(rownames(sdsummary)=="Depletion_t"),"Std. Error"])
      polygon( x=c(1:n_years,n_years:1), y=c(Mat[,1]+Mat[,2],rev(Mat[,1]-Mat[,2])), col=rgb(1,0,0,0.2), border=NA)
      # Catch
      matplot( cbind(catch_t,report$cpred_t[-n_years]), type="l", log="y", ylab="", xlab="year", main="Catch")
      # Exploitation rate
      matplot( cbind(exploit_t[-n_years],report$exploit_t), type="l", log="y", ylab="", xlab="year", main="Exploitation fraction")
    dev.off()
  }

  ######
  # Run without estimating overdispersion
  ######
  Map = list()
  Map[["log_r"]] = factor(NA)
  Map[["log_sigmac"]] = factor(NA)
  Map[["log_sigmam"]] = factor(NA)

  # Compile
  dyn.load( dynlib(Version_statespace) )
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map)

  # Optimize
  Opt = nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1, "eval.max"=1e4, "iter.max"=1e4))
  Opt[["diagnostics"]] = data.frame( "Est"=Opt$par, "final_gradient"=Obj$gr(Opt$par) )
  Report = Obj$report()
  SD = try(sdreport( Obj ))

  # visualize and save results
  if( all(abs(Opt$diagnostics$final_gradient)<0.01) ){
    Results["without",repI,c("k","RE_final_depletion","sigmam")] = c(Report$k, Report$Depletion_t[n_years]/(n_t[n_years]/k), exp(Report$log_sigmam))
    if( !inherits(SD,"try-error")) Plot_Fn( report=Report, sdsummary=summary(SD), dir=DateDir, name=paste0("Rep_",repI,"_Without_overdispersion.png"))
  }

  ######
  # Run with estimating overdispersion
  ######
  Map = list()
  Map[["log_r"]] = factor(NA)
  Map[["log_sigmac"]] = factor(NA)
  Params[["log_sigmam"]] = log(1)

  # Compile
  dyn.load( dynlib(Version_statespace) )
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map)

  # Optimize
  Opt = nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1, "eval.max"=1e4, "iter.max"=1e4))
  Opt[["diagnostics"]] = data.frame( "Est"=Opt$par, "final_gradient"=Obj$gr(Opt$par) )
  Report = Obj$report()
  SD = try(sdreport( Obj ))

  # visualize and save results
  if( all(abs(Opt$diagnostics$final_gradient)<0.01) ){
    Results["with",repI,c("k","RE_final_depletion","sigmam")] = c(Report$k, Report$Depletion_t[n_years]/(n_t[n_years]/k), exp(Report$log_sigmam))
    if( !inherits(SD,"try-error")) Plot_Fn( report=Report, sdsummary=summary(SD), dir=DateDir, name=paste0("Rep_",repI,"_With_overdispersion.png"))
  }
}

#####################
# Compile results
#####################

HistFn = function(x0, x1, trueval, nbreaks=10, ...){
  xlim = range(pretty(range( c(trueval,x0,x1),na.rm=TRUE )))
  Hist0 = hist(x0, plot=FALSE, breaks=seq(xlim[1],xlim[2],length=nbreaks))
  Hist1 = hist(x1, plot=FALSE, breaks=seq(xlim[1],xlim[2],length=nbreaks))
  ylim = c(0, max(c(Hist0$counts,Hist1$counts)) )
  plot(Hist0, xlim=xlim, col=rgb(0,0,0,0.2), ylim=ylim, ...)
  plot(Hist1, xlim=xlim, col=rgb(1,0,0,0.2), ylim=ylim, add=TRUE)
  abline(v=c(mean(x0,na.rm=TRUE),mean(x1,na.rm=TRUE),trueval), lwd=3, col=c("black","red","blue"))
}
png( paste0(DateDir,"_Boxplot_summary.png"), width=10, height=5, res=200, units="in")
  par( mfrow=c(1,3), mar=c(3,3,4,0), mgp=c(2,0.5,0), tck=-0.02, yaxs="i")
  HistFn( x0=Results["with",,"k"], x1=Results["without",,"k"], trueval=k, xlab="", ylab="", main="Carrying capacity", yaxt="n" )
  HistFn( x0=Results["with",,"RE_final_depletion"], x1=Results["without",,"RE_final_depletion"], trueval=1, xlab="", ylab="", main="Error in final\nrelative abundance", yaxt="n" )
  HistFn( x0=Results["with",,"sigmam"], x1=Results["without",,"sigmam"], trueval=sigmam, xlab="", ylab="", main="Overdispersion", yaxt="n" )
dev.off()




