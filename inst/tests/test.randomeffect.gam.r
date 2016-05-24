    

  RLibrary( "mgcv", "chron", "lattice"  ) 
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 


    require(fields)
    require(mgcv)
    require (boot)  # inv.logit
    

  p$vars.to.model = c("R0.mass")
  p$regions = c("cfa4x")
   
  p$model.type = "gam.full" # choose method for habitat model :
  p$habitat.threshold.quantile = 0.05 # quantile at which to consider zero-valued abundance
  p$optimizers = c( "nlm", "perf", "bam", "bfgs", "newton", "Nelder-Mead" )  # used by GAM
	p$prediction.dyear = 0.8 # predict for ~ Sept 1 
  p$threshold.distance = 15  # limit to extrapolation/interpolation in km
     

  v = p$vars.to.model

  set = habitat.model.db( DS="basedata", p=p, v=v )
  xset = set[which(set$yr %in% 2005:2007),]

  ops = c( "outer", "nlm" ) 
  
  # Q0 = gam( Y~s(plon,plat, by=yr), data=xset, optimizer=ops, family=binomial(), select=T)
  
  Q1 = gam( Y~s(plon,plat, by=as.factor(yr), bs="ts" ) + as.factor(yr) + s(t), data=xset, optimizer=ops, family=binomial(), select=T)
  
  Q1a = gam( Y~s(plon,plat, by=as.factor(yr), bs="ts", k=30 ) + as.factor(yr) + s(t), data=xset, optimizer=ops, family=binomial(), select=T)
  Q1b = gam( Y~s(plon,plat, by=as.factor(yr), bs="ts", k=50 ) + as.factor(yr) + s(t), data=xset, optimizer=ops, family=binomial(), select=T)
  Q1c = gam( Y~s(plon,plat, by=as.factor(yr), bs="ts", k=100 ) + as.factor(yr) + s(t), data=xset, optimizer=ops, family=binomial(), select=T)
  Q1d = gam( Y~s(plon,plat, by=as.factor(yr), bs="ts", k=200 ) + as.factor(yr) + s(t), data=xset, optimizer=ops, family=binomial(), select=T)
 
  Q1e = gam( Y~s(plon,plat, by=as.factor(yr), bs="ts", k=256 ) + as.factor(yr) + s(t), data=xset, optimizer=ops, family=binomial(), select=T)
 
  # -- k =256 is too much k=200 seems best (based upon AIC)
 AIC(Q0, Q1, Q1a, Q1b, Q1c, Q1d)
           df      AIC
Q0   62.31907 1323.373  k=internal
Q1   62.31907 1323.373  k=internal
Q1a  62.31907 1323.373  k=50
Q1b  81.06555 1314.660  k=100
Q1c 189.32528 1212.943  k=200 -- sweet spot
Q1d 214.52285 1256.456  k=256 -- not enough data in some years 


# gamm performance is poor  .. possibly due to insufficient data ?
  Q2 = gamm( Y~s(plon,plat, by=as.factor(yr),bs="ts", k=100)+ as.factor(yr) + s(t), correlation=corGaus( c( 75.0, 0.2 ),form=~plon+plat|yr, nugget=TRUE), data=xset, optimizer=ops, family=binomial(), select=T)
  Q2a = gamm( Y~s(plon,plat, by=as.factor(yr),bs="ts" )+ as.factor(yr) , correlation=corGaus( c( 75.0, 0.2 ),form=~plon+plat|yr, nugget=TRUE), data=xset, optimizer=ops, family=binomial(), select=T)
 
  Q2b = gam( Y~s(plon,plat, by=as.factor(yr), bs="re", k=100 ) + as.factor(yr) + s(t), data=xset, optimizer=ops, family=binomial(), select=T)

  AIC(Q0, Q1)



  Q = Q1
  
  summary(Q); 
  AIC (Q) #207506
  plot(Q, all.terms=T, rug=T, jit=T, seWithMean=T, pers=T, trans=inv.logit, scale=0 )
  
  ppp = predict( Q, set, type="response")
  cor(ppp,set$Y, use="pairwise.complete.obs")^2 #.54
  



  # posterior simulations
  Hmodel = habitat.model.db( DS="habitat", v= v )
  Hsim = gam.simulation( M=Hmodel, X= PS, nsims=p$nsims ) #~8min
  rm( Hmodel); gc()
  
  oops = which( is.na(Hsim) ) 
  if (length(oops) > 0)  Hsim[oops ] = 0  # assume to be zero
  #Hsim = round(Hsim)  # force binary


  PS = habitat.db ( DS="complete", year=y, p=p )
  PS$dyear = p$prediction.dyear  # must be same as above
  PS$t = NA
   
  PST = temperature.db( p=p, DS="spatial.interpolation", yr=y  )
  if (is.null(PST)) next ()
 
  dyears = (c(1:(p$nw+1))-1)  / p$nw # intervals of decimal years... fractional year breaks
  dyr = as.numeric( cut( p$prediction.dyear, breaks=dyears, include.lowest=T, ordered_result=TRUE ) ) # integer representation of season

  PS$t = PST[, dyr ]
  PS$t[ which(PS$t < -2) ] = -2
  PS$t[ which(PS$t > 30) ] = 30 

  iitna = which( ! is.finite( PS$t ) ) 
  if (length(iitna)>0) PS$t[iitna] = PS$tmean[iitna]

  PS$z = log(PS$z)
  PS$dt.seasonal = PS$tmean - PS$t 
  PS$dt.annual = PS$tmean - PS$tmean.cl
  PS$sa = 1




  
