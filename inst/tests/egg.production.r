 
  # moved into the main body of the snow crab analytical routines
  # look at:
  #
  # fecundity.indirect()
  # figure.timeseries.fecundity ()
  #
  # both of which are run at:
  #
  # 4.figures.r

   	
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  

  
  # ------------------------------------------
  # krige relevant variable to whole-shelf scale
   
   p$do.parallel = T
   p$vars.to.model =  c(
    "f7.no","m7.no", 
    "f8.no","m8.no", 
    "f9.no","m9.no", 
    "f10.no","m10.no", 
    "totno.female.primiparous", "totno.female.multiparous", 
    "fecundity" 
   )
   p$regions.to.model = "cfaall"

   p$ofname = file.path( project.datadirectory("snowcrab"), "R", "egg.production.rdata" )
   
   p = make.list( list(p$vars.to.model, p$years.to.model), Y=p)


  ###################
  # krige
    krige.ordinary( p, init.files=init.files  )  # ~ 6 hr on "io X 8"
    K = krige.block ( p, init.files=init.files  ) # to est CI and get data summaries ~ 2 days
    krige.map ( p, init.files=init.files  )  # ~ 1 day?
    gmt.cleanup() # clean up any stragglers
    cmd( "rsync -avzn tethys:/home/jae/ecomod/snowcrab/R /home/jae/ecomod/snowcrab/" )



  # ------------------------------------------
  # Fecundity estimated indirectly via total number of females of primi and muli and applying mean egg production from allometric estimate
  fec.ind = fecundity.indirect(  outdir=file.path( p$annual.results, "timeseries", "survey" ) )
  fec.ind$yr = fec.ind$years

 # ------------------------------------------
  # Fecundity estimated directly via kriging and individual-based fecundity estimate
  # p$ofname="/home/jae/ecomod/snowcrab/R/females.kriged.results.rdata"
  
  fec.krig = figure.timeseries.fecundity(p,  outdir=file.path( p$annual.results, "timeseries",  "kriged"  ) )
  fec.krig2 = tapply( fec.krig$total, fec.krig$yr, sum, na.rm=T) 
  fec.krig2l = tapply( fec.krig$lbound, fec.krig$yr, max, na.rm=T)   # this ignores the other areas .. but a good approximation
  fec.krig2u = tapply( fec.krig$ubound, fec.krig$yr, sum, na.rm=T) 

  out.k = as.data.frame(cbind( fec.krig2, fec.krig2l, fec.krig2u ))
  out.k = out.k * 10^6  # return from millions to unit of numbers converted in functions.trawl at add.to.set

  names(out.k)= c("total.egg.kriged", "total.egg.kriged.upperbound", "total.egg.kriged.lowerbound")
  out.k$yr = as.numeric( rownames( out.k) )
  out = merge (fec.ind[, c("yr", "total.egg.all", "total.egg.all.sd")],out.k, by="yr" )  # log10(n)
  x11()
  plot(out$yr, out$total.egg.kriged)


  ############ estimation of natural mortality from egg stage to Instar 8  
  # ... there exist catachability issues with Instar 8 (20-40% undersampled) and so these estimates are conservative

  # Assumption: egg released in year Y are captured as instar 8 crab in year Y+7 (assumes 1 year of egg retention/development by females)
  # Assumption catachbility of Instar 8 = 1
  # Assumption exponential rate of mortality
    

  A = fecundity.indirect(  outdir=file.path( p$annual.results, "timeseries", "survey" ) )
  A$yr = A$years + 7
  A = A[ which(A$years >= 1997) ,]
  #   p$ofname = "/home/jae/ecomod/snowcrab/R/instar8.rdata"
 
  load (p$ofname)
  B = K

  C = merge( B, A, by="yr", all.x=T, all.y=F )
  C = C[ is.finite( C$total + C$total.egg.all ) ,]

  fem = C[ which(C$vars=="f8.no") ,]
  mal = C[ which(C$vars=="m8.no") ,]

  plot( fem$total.egg.all, fem$total)
  lines( mal$total.egg.all, mal$total, pch=21  )

  
  plot( mal$yr, mal$total.egg.all/max(mal$total.egg.all), type="b", ylim=c(0,1) )
  lines ( mal$yr, mal$total/max(mal$total), col="blue" )
  lines( fem$yr, fem$total/max(fem$total), col="red" )
 

  plot( A$years, A$total.egg.all, type="b", ylim=c(0,4e8))
  lines ( mal$yr, mal$total, col="blue" )
  lines( fem$yr, fem$total, col="red" )


 
  plot( A$years, A$total.egg.all/max(A$total.egg.all), type="b", ylim=c(0,1) )
#   lines ( mal$yr, mal$total/max(mal$total), col="blue" )
  lines( fem$yr, fem$total/max(fem$total), col="red" )
 
 hh = fem[ , c("yr", "total", "lbound", "ubound")]
 ii = A[,c("years", "total.egg.all","total.egg.all.sd", "upper", "lower" )]

 out = merge(hh, ii, by.x="yr", by.y="years")
 out$total.egg.all.lower = 10^out$lower
 out$total.egg.all.upper = 10^out$upper
 out$upper = NULL
 out$lower = NULL

 out$instar8 = out$total
 out$total = NULL
 out$lbound = NULL
 out$ubound = NULL



# mortality rates estimates crude:
 lep = 190000; rt=4; mt=8.5 
 Z=(log( lep/2 ) - log(1) )/ (rt + mt - 1)
 S = exp(-Z)
 A = 1-S
[1] 0.369109

# population/survey-based:
Z = (log(fem$total.egg.all/2 ) - log(fem$total)) / 7
S = exp( - Z)
A = 1 - S

mean(A)
exp( mean(A))



yr
   1997    1998    1999    2000 
1.70772 1.68391 1.60439 1.59819 

> mean ((log(fem$total.egg.all/2 ) - log(fem$total)) / 7)
[1]  1.64855

> sd ((log(fem$total.egg.all/2 ) - log(fem$total)) / 7)
[1] 0.0554936


# variance estimate from sets :: kriging estimates are unstable

s=snowcrab.db("set")
x = log( s$f8.no + exp(1) )

f8 = NULL
f8$n = length(x)
f8$mean = mean( x )
f8$sd =  sd( x)
f8$se = f8$sd/ sqrt(f8$n-1)

u = rlnorm( 1000, f8$mean, f8$sd )

out = NULL
for (i in 1: nrows(fem) {
  


}


