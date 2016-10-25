
############################################################################
# Compile model cpp file in TMB
############################################################################
library(TMB);
compile("biomassdynamic.cpp"); 
dyn.load(dynlib("biomassdynamic"));

###########################################################################
# Estimation
###########################################################################
# The objective function
obj <- MakeADFun(data, parameters,                 # <<< data and initial values
                 random=c("log_P","log_r","log_K"),# <<< specify random variables
                 DLL="biomassdynamic")                    # <<< compiled model file
# Optimize the objective function
opt <- nlminb(obj$par,obj$fn,obj$gr)
# Report estimates and standard errors
rep <- sdreport(obj)

