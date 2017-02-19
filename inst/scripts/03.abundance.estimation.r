

# --------------------------------------------------------------
#  Ensure the following scripts complete without error:
#  these external dependencies permit lookup of data, for snowcrab.db("complete.redo") :
#  system.file(package="bio.indicators", "scripts", "01.indicators.r")  
#  system.file(package="bio.indicators", "scripts", "02.*.r")  


# 1. Define some additional starting parameters for debugging
#    choose various over-rides: these are initially defined in parameters.r

current.year = 2016
p = bio.snowcrab::load.environment( year.assessment=current.year )




# --------------------------------------------------------------
# using environmental data ... estimate/lookup missing environmental data .. (t,z)
logbook.db( DS ="fisheries.complete.redo", p=p )
snowcrab.db( DS ="set.complete.redo", p=p )



# -------------------------------------------------------------------------------------
# abundance .. positive valued data .. vn = "snowcrab.large.males_abundance"

p = bio.snowcrab::load.environment( year.assessment=current.year )
p$selection=list( 
  name = "snowcrab.large.males_abundance",
  type = "abundance",
  sex=0, # male
  mat=1, # do not use maturity status in groundfish data as it is suspect .. 
  spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
  len= c( 95, 200 )/10, #  mm -> cm ; indicators.db in cm
  drop.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable 
)
p$lbm_local_modelengine = "twostep"
# p$lbm_global_family = gaussian(link="log")
# p$lbm_local_family = gaussian(link="log") 

# 11 hrs with these settings
p$lbm_twostep_space = "krige"
p$lbm_gam_optimizer=c("outer", "bfgs") 
p$lbm_distance_statsgrid = 2 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
p$lbm_distance_prediction = p$lbm_distance_statsgrid * 0.75 # this is a half window km
p$lbm_distance_scale = 50

p = bio.snowcrab::snowcrab.parameters( p=p, DS="lbm", varname=p$selection$name  )

# o = snowcrab_lbm(p=p, DS="lbm_inputs" )  # create fields for 
DATA='snowcrab_lbm( p=p, DS="lbm_inputs" )'
lbm( p=p, DATA=DATA, tasks=c("initiate", "globalmodel") ) # 30 min
#   lbm( p=p, tasks=c( "stage0" ) ) # serial mode
#   lbm( p=p, tasks=c( "continue" ) )    
lbm( p=p, tasks=c( "stage1" ) ) #  3 hrs 
lbm( p=p, tasks=c( "stage2" ) ) #   1 hrs
lbm( p=p, tasks=c( "save" ) )

p = make.list( list( yrs=p$yrs), Y=p )
parallel.run( snowcrab_lbm, p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_lbm( p=p, DS="lbm.stats.redo" ) # warp stats to other grids
snowcrab_lbm( p=p, DS="complete.redo" )
snowcrab_lbm( p=p, DS="baseline.redo" )
snowcrab_lbm( p=p, DS="map.all" )

global_model = lbm_db( p=p, DS="global_model") 
summary( global_model )
plot(global_model)


Family: gaussian 
Link function: identity 

Formula:
snowcrab.large.males_abundance ~ s(t, k = 3, bs = "ts") + s(tmean.climatology, 
    k = 3, bs = "ts") + s(tsd.climatology, k = 3, bs = "ts") + 
    s(log(z), k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + 
    s(log(ddZ), k = 3, bs = "ts") + s(log(mr), k = 3, bs = "ts") + 
    s(Npred, k = 3, bs = "ts") + s(smr, k = 3, bs = "ts") + s(log.substrate.grainsize, 
    k = 3, bs = "ts") + s(ca1, k = 3, bs = "ts") + s(ca2, k = 3, 
    bs = "ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  6.70334    0.01501   446.5   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                              edf Ref.df       F  p-value    
s(t)                       1.0690      2   9.527 4.53e-06 ***
s(tmean.climatology)       1.5693      2   4.250  0.00518 ** 
s(tsd.climatology)         1.8496      2  10.253 1.58e-05 ***
s(log(z))                  2.0000      2 125.700  < 2e-16 ***
s(log(dZ))                 0.5762      2   0.911  0.05426 .  
s(log(ddZ))                0.4980      2   0.725  0.06517 .  
s(log(mr))                 1.9061      2 183.476  < 2e-16 ***
s(Npred)                   0.8951      2   5.169  0.00064 ***
s(smr)                     1.9750      2  31.648 1.24e-14 ***
s(log.substrate.grainsize) 1.9324      2  72.260  < 2e-16 ***
s(ca1)                     2.0000      2  26.372 1.57e-12 ***
s(ca2)                     1.5737      2 148.871  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.216   Deviance explained = 21.8%
GCV = 1.0687  Scale est. = 1.0644    n = 4722


# -------------------------------------------------
# presence-absence
# current.year = 2016

p = bio.snowcrab::load.environment( year.assessment=current.year )
p$selection=list( 
  name = "snowcrab.large.males_presence_absence",
  type = "presence_absence",
  sex=0, # male
  mat=1, # do not use maturity status in groundfish data as it is suspect ..   
  spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
  len= c( 95, 200 )/10, #  mm -> cm ; indicators.db in cm
  drop.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable 
)

p$lbm_global_family = binomial()

p$lbm_local_modelengine = "twostep"
p$lbm_local_family = gaussian()  # after logit transform by global model, it becomes gaussian (logit scale)
p$lbm_twostep_space = "krige"
p$lbm_gam_optimizer=c("outer", "bfgs") 
p$lbm_distance_statsgrid = 2 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
p$lbm_distance_prediction = p$lbm_distance_statsgrid * 0.75 # this is a half window km
p$lbm_distance_scale = 50


p = bio.snowcrab::snowcrab.parameters( p=p, DS="lbm", varname=p$selection$name  )

# o = snowcrab_lbm(p=p, DS="lbm_inputs" )  # create fields for 
DATA='snowcrab_lbm( p=p, DS="lbm_inputs" )'

lbm( p=p, DATA=DATA, tasks=c("initiate", "globalmodel") ) # 30 min
#   lbm( p=p, tasks=c( "stage0" ) ) # serial mode
#   lbm( p=p, tasks=c( "continue" ) )    
lbm( p=p, tasks=c( "stage1" ) ) #  3 hrs 
lbm( p=p, tasks=c( "stage2" ) ) #   1 hrs
lbm( p=p, tasks=c( "save" ) )

p = make.list( list( yrs=p$yrs), Y=p )
parallel.run( snowcrab_lbm, p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_lbm( p=p, DS="lbm.stats.redo" ) # warp stats to other grids
snowcrab_lbm( p=p, DS="complete.redo" )
snowcrab_lbm( p=p, DS="baseline.redo" )
snowcrab_lbm( p=p, DS="map.all" )

global_model = lbm_db( p=p, DS="global_model") 
summary( global_model )
plot(global_model)

Family: binomial 
Link function: logit 

Formula:
snowcrab.large.males_presence_absence ~ s(t, k = 3, bs = "ts") + 
    s(tmean.climatology, k = 3, bs = "ts") + s(tsd.climatology, 
    k = 3, bs = "ts") + s(log(z), k = 3, bs = "ts") + s(log(dZ), 
    k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts") + s(log(mr), 
    k = 3, bs = "ts") + s(Npred, k = 3, bs = "ts") + s(smr, k = 3, 
    bs = "ts") + s(log.substrate.grainsize, k = 3, bs = "ts") + 
    s(ca1, k = 3, bs = "ts") + s(ca2, k = 3, bs = "ts")

Parametric coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  0.82108    0.03282   25.02   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                                 edf Ref.df  Chi.sq  p-value    
s(t)                       1.769e+00      2 105.703  < 2e-16 ***
s(tmean.climatology)       1.838e+00      2  68.113  < 2e-16 ***
s(tsd.climatology)         1.910e+00      2  14.130 0.000585 ***
s(log(z))                  1.569e+00      2 170.503  < 2e-16 ***
s(log(dZ))                 6.303e-01      2   2.517 0.035988 *  
s(log(ddZ))                9.639e-01      2  12.371 6.30e-05 ***
s(log(mr))                 1.855e+00      2 196.086  < 2e-16 ***
s(Npred)                   7.202e-05      2   0.000 0.732550    
s(smr)                     9.897e-01      2  15.089 3.94e-05 ***
s(log.substrate.grainsize) 1.752e+00      2 213.989  < 2e-16 ***
s(ca1)                     1.849e+00      2  22.419 4.75e-06 ***
s(ca2)                     1.857e+00      2 246.368  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.296   Deviance explained = 25.4%
UBRE = -0.031615  Scale est. = 1         n = 6853






# collect all predictions into a single file and return:
p = bio.snowcrab::load.environment( year.assessment=current.year )
p$selection=list( 
  name = "snowcrab.large.males_abundance",
  type = "abundance",
  sex=0, # male
  mat=1, # do not use maturity status in groundfish data as it is suspect .. 
  spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
  len= c( 95, 200 )/10, #  mm -> cm ; indicators.db in cm
  drop.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable 
}
p = bio.snowcrab::snowcrab.parameters( p=p, DS="lbm", varname=p$selection$name  )

interpolation.db( DS="biomass.redo", p=p  )
interpolation.db( DS="biomass.map", p=p  )

K = interpolation.db( DS="timeseries", p=p  )
str(K)

# table.view( K )
# figure.timeseries.errorbars( K[], outdir=outdir, fname=paste(vv, rr, sep=".") )


# update data summaries of the above results
p$vars.tomodel="R0.mass"
biomass.summary.db("complete.redo", p=p) #Uses the model results to create a habitat area expanded survey index

# biomass.summary.survey.db("complete.redo", p=p)#Uses average surface area from the past 5 years if a habitat area expanded surface area is not possible -- JC .. no longer used .. marked for deletion




### --------- prediction success:
set = snowcrab_lbm(p=p, DS="input_data", voi=p$selection$name )

S = set[ , c("plon", "plat") ]

ii = array_map( "xy->1", S, gridparams=p$gridparams )
bs = bathymetry.db(p=p, DS="baseline")
bb = array_map( "xy->1", bs, gridparams=p$gridparams )
im = match(  ii, bb )
it = match( set$yr, p$yrs )

bm = interpolation.db( DS="biomass", p=p  )
spred = bm$m[cbind(im, it)]  # approximate match (ignoring seasonality)

summary ( lm(spred~snowcrab.large.males_abundance, data=set, na.actio="na.omit" ) )
plot(spred~snowcrab.large.males_abundance, data=set )
cor(spred,set$snowcrab.large.males_abundance, use="complete.obs")


# determine presence absence(Y) and weighting(wt)
#      set$weekno = floor(set$julian / 365 * 52) + 1
#      set$dyear = floor(set$julian / 365 ) + 1


