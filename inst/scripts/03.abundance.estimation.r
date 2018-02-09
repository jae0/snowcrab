
  if (!exists("year.assessment")) {
    year.assessment=lubridate::year(Sys.Date()) -1
    year.assessment=lubridate::year(Sys.Date())
  }

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )



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

    p = stmv( p=p, runmode=c("initialize", "globalmodel" ), use_saved_state=FALSE ) # no global_model and force a clean restart

    currentstatus = stmv_db( p=p, DS="statistics.status" )
    p = parallel_run( stmv_interpolate, p=p, 
      runindex=list( locs=currentstatus$todo[sample.int(length( currentstatus$todo ))] ),
      local.n.complete=currentstatus["n.complete"]
    ) 
    stmv_db( p=p, DS="save_current_state" ) # saved current state (internal format)
    if (exists("cl", p)) stopCluster( p$cl )
  
    currentstatus = stmv_db(p=p, DS="statistics.status.reset" )
    parallel_run( stmv_interpolate, p=p, 
      runindex=list( locs=currentstatus$todo[sample.int(length( currentstatus$todo ))] ),
      local.n.complete=currentstatus["n.complete"], 
      stmv_distance_max=p$stmv_distance_max*mult, 
      stmv_distance_scale=p$stmv_distance_scale*mult 
    )
    stmv_db( p=p, DS="save_current_state" ) # saved current state 
    if (exists("cl", p)) stopCluster( p$cl )

    # currentstatus = stmv_db( p=p, DS="statistics.status.reset" )
    # p = parallel_run( stmv_interpolate, p=p, 
    #   runindex=list( locs= currentstatus$todo[sample.int(length( currentstatus$todo ))] ),
    #   local.n.complete=currentstatus["n.complete"], 
    #   stmv_local_modelengine = "tps" ) 
    # stmv_db( p=p, DS="save_current_state" )
    # if (exists("cl", p)) stopCluster( p$cl )

    stmv_db( p=p, DS="stmv.prediction.redo" ) # save to disk for use outside stmv*, returning to user scale
    stmv_db( p=p, DS="stats.to.prediction.grid.redo") # save to disk for use outside stmv*
  
    # if (really.finished) stmv_db( p=p, DS="cleanup.all" )



snowcrab_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids
snowcrab_stmv( p=p, DS="complete.redo" )
snowcrab_stmv( p=p, DS="baseline.redo" )
snowcrab_stmv( p=p, DS="map.all" )

global_model = stmv_db( p=p, DS="global_model")
summary( global_model )
plot(global_model)

#Below lines are just model outputs for comparison sake
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
(Intercept)  5.58131    0.01869   298.6   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                                 edf Ref.df       F  p-value
s(t)                       1.461e+00      2  67.664  < 2e-16 ***
s(tmean.climatology)       1.441e+00      2  19.462 6.04e-11 ***
s(tsd.climatology)         1.951e+00      2  15.470 1.17e-07 ***
s(log(z))                  1.548e+00      2 199.031  < 2e-16 ***
s(log(dZ))                 8.326e-01      2   3.413  0.00213 **
s(log(ddZ))                9.484e-01      2   5.210  0.00026 ***
s(log(mr))                 1.507e+00      2 219.813  < 2e-16 ***
s(Npred)                   1.360e-08      2   0.000  0.60121
s(smr)                     1.734e+00      2  17.919 1.88e-09 ***
s(log.substrate.grainsize) 1.898e+00      2 189.463  < 2e-16 ***
s(ca1)                     1.947e+00      2  48.772  < 2e-16 ***
s(ca2)                     1.901e+00      2 292.569  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.357   Deviance explained = 35.9%
GCV = 2.4006  Scale est. = 2.3943    n = 6853

# -------------------------------------------------
# presence-absence
# year.assessment = 2016

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

    p = stmv( p=p, runmode=c("initialize", "globalmodel" ), use_saved_state=FALSE ) # no global_model and force a clean restart

    currentstatus = stmv_db( p=p, DS="statistics.status" )
    p = parallel_run( stmv_interpolate, p=p, 
      runindex=list( locs=currentstatus$todo[sample.int(length( currentstatus$todo ))] ),
      local.n.complete=currentstatus["n.complete"]
    ) 
    stmv_db( p=p, DS="save_current_state" ) # saved current state (internal format)
    if (exists("cl", p)) stopCluster( p$cl )
  
    currentstatus = stmv_db(p=p, DS="statistics.status.reset" )
    parallel_run( stmv_interpolate, p=p, 
      runindex=list( locs=currentstatus$todo[sample.int(length( currentstatus$todo ))] ),
      local.n.complete=currentstatus["n.complete"], 
      stmv_distance_max=p$stmv_distance_max*mult, 
      stmv_distance_scale=p$stmv_distance_scale*mult 
    )
    stmv_db( p=p, DS="save_current_state" ) # saved current state 
    if (exists("cl", p)) stopCluster( p$cl )

    # currentstatus = stmv_db( p=p, DS="statistics.status.reset" )
    # p = parallel_run( stmv_interpolate, p=p, 
    #   runindex=list( locs= currentstatus$todo[sample.int(length( currentstatus$todo ))] ),
    #   local.n.complete=currentstatus["n.complete"], 
    #   stmv_local_modelengine = "tps" 
    # ) 
    # stmv_db( p=p, DS="save_current_state" )
    # if (exists("cl", p)) stopCluster( p$cl )

    stmv_db( p=p, DS="stmv.prediction.redo" ) # save to disk for use outside stmv*, returning to user scale
    stmv_db( p=p, DS="stats.to.prediction.grid.redo") # save to disk for use outside stmv*
  
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
