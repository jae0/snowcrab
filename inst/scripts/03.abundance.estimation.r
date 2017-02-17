

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
p$lbm_distance_statsgrid = 3 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
p$lbm_distance_prediction = p$lbm_distance_statsgrid*0.75  # this is a half window km
p$lbm_distance_scale = 25

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
    s(log(dZ), k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts") + 
    s(log(mr), k = 3, bs = "ts") + s(Npred, k = 3, bs = "ts") + 
    s(smr, k = 3, bs = "ts") + s(log.substrate.grainsize, k = 3, 
    bs = "ts") + s(ca1, k = 3, bs = "ts") + s(ca2, k = 3, bs = "ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)   5.5813     0.0192   290.6   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                                 edf Ref.df      F  p-value    
s(t)                       1.855e+00      2 107.71  < 2e-16 ***
s(tmean.climatology)       1.665e+00      2  31.49  < 2e-16 ***
s(tsd.climatology)         1.775e+00      2  14.23 1.35e-07 ***
s(log(dZ))                 1.746e+00      2  10.52 6.34e-06 ***
s(log(ddZ))                1.637e+00      2  11.79 1.04e-06 ***
s(log(mr))                 1.727e+00      2 229.51  < 2e-16 ***
s(Npred)                   1.212e-08      2   0.00    0.751    
s(smr)                     1.499e+00      2  19.02 2.19e-10 ***
s(log.substrate.grainsize) 1.847e+00      2 287.68  < 2e-16 ***
s(ca1)                     1.916e+00      2  37.58  < 2e-16 ***
s(ca2)                     1.946e+00      2 110.24  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.321   Deviance explained = 32.3%
GCV = 2.5344  Scale est. = 2.5275    n = 6853



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
p$lbm_distance_statsgrid = 3 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
p$lbm_distance_prediction = p$lbm_distance_statsgrid*0.75 # this is a half window km
p$lbm_distance_scale = 25


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
    k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + s(log(ddZ), 
    k = 3, bs = "ts") + s(log(mr), k = 3, bs = "ts") + s(Npred, 
    k = 3, bs = "ts") + s(smr, k = 3, bs = "ts") + s(log.substrate.grainsize, 
    k = 3, bs = "ts") + s(ca1, k = 3, bs = "ts") + s(ca2, k = 3, 
    bs = "ts")

Parametric coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  0.78170    0.03176   24.61   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                             edf Ref.df  Chi.sq  p-value    
s(t)                       1.937      2 159.881  < 2e-16 ***
s(tmean.climatology)       1.887      2 109.616  < 2e-16 ***
s(tsd.climatology)         1.680      2  10.115  0.00270 ** 
s(log(dZ))                 1.255      2   9.449  0.00177 ** 
s(log(ddZ))                1.512      2  20.517 3.52e-06 ***
s(log(mr))                 1.342      2 158.152  < 2e-16 ***
s(Npred)                   1.903      2  20.863 1.58e-05 ***
s(smr)                     1.737      2   8.086  0.00968 ** 
s(log.substrate.grainsize) 1.597      2 298.690  < 2e-16 ***
s(ca1)                     1.727      2  20.945 6.30e-06 ***
s(ca2)                     1.918      2 116.735  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.271   Deviance explained =   23%
UBRE = -0.001273  Scale est. = 1         n = 6853
---



# collect all predictions into a single file and return:
interpolation.db( DS="biomass.redo", p=p  )
interpolation.db( DS="biomass.map", p=p  )

K = interpolation.db( DS="timeseries", p=p  )
str(K)

# table.view( K )
# figure.timeseries.errorbars( K[], outdir=outdir, fname=paste(vv, rr, sep=".") )


# update data summaries of the above results
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


