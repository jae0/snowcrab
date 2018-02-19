
  if (!exists("year.assessment")) {
    year.assessment=lubridate::year(Sys.Date()) -1
    year.assessment=lubridate::year(Sys.Date())
  }

# --------------------------------------------------------------
#  Ensure the following scripts complete without error:
#  these external dependencies permit lookup of data, for this script
#Require the following:
#snowcrab.db("complete.redo")

#Run aegis::(inst/scripts/11.surveys.r), need to update year within or run next line
  #system.file(package="aegis", "scripts", "11.surveys.r")

#Run aegis::(inst/scripts/temperature.r), need to update year within
  #system.file(package="aegis", "scripts", "07.temperature.R") or run next line

#Substrate and bathymetry can be run (as below) if suspect significant changes in one or both of these datasets

#BZ 2017 these lines below can directly run the indicators without goint to run individual scripts
#  system.file(package="aegis", "scripts", "07.temperature.R")
#  system.file(package="aegis", "scripts", "10.lookuptables.r")
#  system.file(package="aegis", "scripts", "11.surveys.r")
#  system.file(package="aegis", "scripts", "16.speciescomposition.R")

# 1. Define some additional starting parameters for debugging
#    choose various over-rides: these are initially defined in parameters.r


# --------------------------------------------------------------
# using environmental data ... estimate/lookup missing environmental data .. (t,z)
#BZ 2017 below lines shouldn't be required. datasets created in 01.snowcrab
#logbook.db( DS ="fisheries.complete.redo", p=p )
#snowcrab.db( DS ="set.complete.redo", p=p )

# -------------------------------------------------------------------------------------
# abundance .. positive valued data .. vn = "snowcrab.large.males_abundance"
# year.assessment = 2017
p = bio.snowcrab::load.environment( year.assessment=year.assessment ) 

# 11 hrs with these settings,
p = snowcrab_stmv( p=p, DS="parameters",
  selection=list(
    name = "snowcrab.large.males_abundance",
    type = "abundance",
    sex=0, # male
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
    len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
    drop.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
  ),
  DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs" )',
  stmv_global_family = gaussian(link=log),
  stmv_local_modelengine = "twostep",
  stmv_twostep_space = "krige",
  stmv_gam_optimizer=c("outer", "bfgs") ,
  stmv_distance_statsgrid = 2, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_scale = 50
)

# o = snowcrab_stmv(p=p, DS="stmv_inputs" )  # create fields for  

stmv( p=p, runmode=c("globalmodel", "interpolate" )  ) # no global_model and force a clean restart

# if (really.finished) stmv_db( p=p, DS="cleanup.all" )

snowcrab_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids (if any)
snowcrab_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids (if any)
snowcrab_stmv( p=p, DS="complete.redo" )
snowcrab_stmv( p=p, DS="baseline.redo" )
snowcrab_stmv( p=p, DS="map.all" )

global_model = stmv_db( p=p, DS="global_model")
summary( global_model )
plot(global_model)

#Below lines are just model outputs for comparison sake
Family: gaussian 
Link function: log 

Formula:
snowcrab.large.males_abundance ~ s(t, k = 3, bs = "ts") + s(tmean.climatology, 
    k = 3, bs = "ts") + s(tsd.climatology, k = 3, bs = "ts") + 
    s(log(z), k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + 
    s(log(ddZ), k = 3, bs = "ts") + s(log.substrate.grainsize, 
    k = 3, bs = "ts") + s(pca1, k = 3, bs = "ts") + s(pca2, k = 3, 
    bs = "ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.64886    0.00403     409   <2e-16

Approximate significance of smooth terms:
                            edf Ref.df      F p-value
s(t)                       1.97      2  39.02 < 2e-16
s(tmean.climatology)       1.64      2   2.78 0.03376
s(tsd.climatology)         1.97      2 191.60 < 2e-16
s(log(z))                  1.54      2 146.66 < 2e-16
s(log(dZ))                 1.99      2   7.77 0.00037
s(log(ddZ))                1.08      2  13.14 5.9e-08
s(log.substrate.grainsize) 2.00      2  56.62 < 2e-16
s(pca1)                    2.00      2 147.26 < 2e-16
s(pca2)                    1.96      2 135.94 < 2e-16

R-sq.(adj) =  0.365   Deviance explained = 36.6%
GCV = 0.011049  Scale est. = 0.011023  n = 7255
 

# -------------------------------------------------
# presence-absence
# year.assessment = 2017
p = bio.snowcrab::load.environment( year.assessment=year.assessment )
p = snowcrab_stmv( p=p, DS="parameters",
  selection=list(
    name = "snowcrab.large.males_presence_absence",
    type = "presence_absence",
    sex=0, # male
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
    len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
    drop.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
  ),
  DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs" )',
  stmv_global_family = binomial(),
  stmv_local_modelengine = "twostep",
  stmv_twostep_space = "krige",
  stmv_gam_optimizer=c("outer", "bfgs"),
  stmv_distance_statsgrid = 2, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** ),
  stmv_distance_scale = 50
)
# o = snowcrab_stmv(p=p, DS="stmv_inputs" )  # create fields for
stmv( p=p, runmode=c("globalmodel", "interpolate" ) ) # no global_model and force a clean restart

# stmv_db( p=p, DS="stmv.results" ) # save to disk for use outside stmv*, returning to user scale
# if (really.finished) stmv_db( p=p, DS="cleanup.all" )

snowcrab_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids
snowcrab_stmv( p=p, DS="complete.redo" )
snowcrab_stmv( p=p, DS="baseline.redo" )
snowcrab_stmv( p=p, DS="map.all" )

global_model = stmv_db( p=p, DS="global_model")
summary( global_model )
plot(global_model)
Family: binomial 
Link function: logit 

Formula:
snowcrab.large.males_presence_absence ~ s(t, k = 3, bs = "ts") + 
    s(tmean.climatology, k = 3, bs = "ts") + s(tsd.climatology, 
    k = 3, bs = "ts") + s(log(z), k = 3, bs = "ts") + s(log(dZ), 
    k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts") + s(log.substrate.grainsize, 
    k = 3, bs = "ts") + s(pca1, k = 3, bs = "ts") + s(pca2, k = 3, 
    bs = "ts")

Parametric coefficients:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)    1.325      0.019    69.6   <2e-16

Approximate significance of smooth terms:
                            edf Ref.df  Chi.sq p-value
s(t)                       2.00      2  410.59  <2e-16
s(tmean.climatology)       1.98      2  181.92  <2e-16
s(tsd.climatology)         2.00      2  436.52  <2e-16
s(log(z))                  2.00      2 2891.37  <2e-16
s(log(dZ))                 1.96      2    1.94   0.370
s(log(ddZ))                1.85      2    9.84   0.005
s(log.substrate.grainsize) 1.98      2  143.32  <2e-16
s(pca1)                    2.00      2  512.02  <2e-16
s(pca2)                    2.00      2  826.64  <2e-16

R-sq.(adj) =  0.563   Deviance explained = 50.3%
UBRE = -0.61879  Scale est. = 1         n = 64488
 


# collect all predictions into a single file and return:
p = bio.snowcrab::load.environment( year.assessment=year.assessment )

p$selection=list(
  name = "snowcrab.large.males_abundance",
  type = "abundance",
  sex=0, # male
  mat=1, # do not use maturity status in groundfish data as it is suspect ..
  spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
  len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
  drop.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
)

p = snowcrab_stmv( p=p, DS="parameters" )

interpolation.db( DS="biomass.redo", p=p  )
interpolation.db( DS="biomass.map", p=p  )

K = interpolation.db( DS="timeseries", p=p  )
str(K)

# table.view( K )

figure.timeseries.snowcrab.habitat(p=p) # /bio.data/bio.snowcrab/assessments/2016/timeseries/interpolated/snowcrab.habitat.sa.png

figure.timeseries.snowcrab.habitat.temperatures(p=p) # /bio.data/bio.snowcrab/assessments/2016/timeseries/interpolated/mean.bottom.temp.snowcrab.habitat.png



# update data summaries of the above results
p$vars.tomodel="R0.mass"
biomass.summary.db("complete.redo", p=p) #Uses the model results to create a habitat area expanded survey index



### --------- prediction success:
set = snowcrab_stmv(p=p, DS="input_data", voi=p$selection$name )

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
