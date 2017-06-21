
require(bio.base)
  
if (!exists("current.year")) current.year=lubridate::year(Sys.Date())
p = bio.snowcrab::load.environment( year.assessment=current.year)

debug.region="cfa4x"
debug.region="cfanorth" 
debug.region="cfasouth"

sb = surplusproduction.db( DS="LaplacesDemon.debug", sourcedata="nosa", debug.region=debug.region) 
#sb = surplusproduction.db( DS="jags.2016", sourcedata="nosa" ) 


# Data
data=list(N=sb$Ndata, IOA=sb$O, CAT=sb$removals);
# Parameters: initial values
parameters=list(log_sigmap=log(0.2),
                log_sigmao=log(0.2),
                log_Q=log(0.1),
                log_r=log(0.8),
                log_K=log(70),
                log_P=rep(log(0.6),sb$Ndata))

#parameters=list(log_tau=rep(log(0.1),sb$U),
#                log_sigma=rep(log(0.1),sb$U),
#                log_Q=rep(log(0.1),sb$U),
#                log_r=rep(log(0.8),sb$U),
#                log_K=log(c(6,70,2)),
#                log_P=matrix(rep(log(0.6),sb$N*sb$U),sb$N,sb$U))
#

############################################################################
# Compile model cpp file in TMB
############################################################################

tmb.dir = file.path(project.codedirectory('bio.snowcrab'),"inst","tmb")
setwd(tmb.dir)

library(TMB);
compile("biomassdynamic.cpp")
dyn.load(dynlib("biomassdynamic"))

###########################################################################
# Estimation
###########################################################################

# The objective function
obj <- MakeADFun(data, parameters,  random="log_P",	DLL="biomassdynamic")  

# Optimize the objective function
opt <- nlminb(obj$par,obj$fn,obj$gr)

# Report estimates and standard errors
rep <- sdreport(obj)

MLE <- list(year=1:data$N,B=rep$value[1:data$N], B.sd=rep$sd[1:data$N],predI=rep$value[1:data$N+data$N], predI.sd=rep$sd[1:data$N+data$N], K=rep$value[data$N*2+1], K.sd=rep$sd[data$N*2+1], r=rep$value[data$N*2+2],r.sd=rep$sd[data$N*2+2],Q=rep$value[data$N*2+3],Q.sd=rep$sd[data$N*2+3],sigma=rep$value[data$N*2+4],sigma.sd=rep$sd[data$N*2+4],tau=rep$value[data$N*2+5],tau.sd=rep$sd[data$N*2+5],nll=rep$value[data$N*2+6],nll.sd=rep$sd[data$N*2+6])
  
  
###########################################################################
# Results plots
###########################################################################
plot(data$IOA,ylim=c(0,max(MLE$predI+MLE$predI.sd)),pch=16,col='red')
with(MLE,lines(year,predI))
with(MLE,lines(year,predI+predI.sd,lty=2))
with(MLE,lines(year,predI-predI.sd,lty=2))


