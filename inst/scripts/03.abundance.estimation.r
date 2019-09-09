require(aegis)

#Pick whichever year reference below is correct (most often year.assessment...-1)
  if (!exists("year.assessment")) {
    year.assessment=lubridate::year(Sys.Date()) -1
    year.assessment=lubridate::year(Sys.Date())
  }

# --------------------------------------------------------------
#  Ensure the following scripts complete without error:
#  these external dependencies permit lookup of data, for this script

#Require the following:
#snowcrab.db("complete.redo")

#BZ Feb 2019- Indented lines can likely be removed as they are repeated with additoinal details below.
            #Run aegis::(inst/scripts/10.surveys.r), need to update year within or run next line
              #system.file(package="aegis", "scripts", "10.surveys.r")

            #Run aegis::(inst/scripts/05.temperature.r), need to update year within
              #system.file(package="aegis", "scripts", "05.temperature.R") or run next line

            #Substrate and bathymetry can be run (as below) if suspect significant changes in one or both of these datasets


#BZ- Jan 2019 Run the following steps before moving to stmv abundance estimation step
# 01.#  system.file(package="aegis", "scripts", "05.temperature.R") #BC Jan 2019- 50 hours with following setings:
    # stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    # stmv_distance_scale = 25, # km ... approx guess of 95% AC range
    # stmv_distance_max = 25*1.5,
    # sampling = c( 1, 1.1, 1.25 ), # fractions of distance scale and n.min to try when insufficient data


#Deprectated 2019. now built into next step (10. surveys) system.file(package="aegis", "scripts", "20.lookuptables.r")


# 02.#  system.file(package="aegis", "scripts", "10.surveys.r")


# 03.#  system.file(package="aegis", "scripts", "11.speciescomposition.R") #February 2018 BZ-~10 hours with following settings:
#drastically reduced required time moving local spatial model to "fft" from GAM
#Needed to run PCA1, reset R, run PCA2, otherwise had a memory issue, didn't dump big.memory
#stmv_rsquared_threshold = 0.2, # lower threshold
#stmv_distance_statsgrid = 4, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
#stmv_distance_scale = c(50, 60, 80), # km ... approx guess of 95% AC range .. data tends to be sprse realtive to pure space models
#  stmv_distance_prediction_max = 4 * 1.25 , # upper limit in distnace to predict upon (just over the grid size of statsgrid) .. in timeseries can become very slow so try to be small


# 1. Define some additional starting parameters for debugging
#    choose various over-rides: these are initially defined in parameters.r

p = bio.snowcrab::load.environment( year.assessment=year.assessment )

# --------------------------------------------------------------
# using environmental data ... estimate/lookup missing environmental data .. (t,z)
#BZ 2017 below lines shouldn't be required. datasets created in 01.snowcrab
#logbook.db( DS ="fisheries.complete.redo", p=p )
#snowcrab.db( DS ="set.complete.redo", p=p )




# -------------------------------------------------------------------------------------
# STEP ONE commercial abundance
#--------------------------------------------------------------------------------------

# abundance .. positive valued data ..
# takes about 5 hrs .. ~1 GB / process
# vn = "snowcrab.large.males_abundance"
# year.assessment = 2017

#----------------------------------------
#Setting up parameters (model formulas, etc)
#----------------------------------------

# 11 hrs with these settings
ncpus = parallel::detectCores()
if (0) {
  ram_required_main_process = ? # GB
  ram_required_per_process  = ?  # about 1.2GB on average ..in 2018, for twostep / fft
  ncpu = min( parallel::detectCores(), trunc( (ram_local()-ram_required_main_process) / ram_required_per_process ) )
}


p = snowcrab_stmv( p=p, DS="parameters",
  variables=list(Y="snowcrab.large.males_abundance"),
  selection=list(
    type = "biomass",
    biologicals=list(
      sex=0, # male
      mat=1, # do not use maturity status in groundfish data as it is suspect ..
      spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
      len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
      ranged_data="len"
    ),
    survey=list(
      data.source = c("snowcrab"),
      yr = p$yrs      # time frame for comparison specified above
    )
  ),
  DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs" )',
  stmv_global_modelengine ="gam",
  stmv_global_family = gaussian(link="log"),
  stmv_local_modelengine = "twostep",
  stmv_twostep_time = "gam",
  stmv_twostep_space = "fft",
  stmv_fft_filter="matern",  #  matern, krige (very slow), lowpass, lowpass_matern
  stmv_gam_optimizer=c("outer", "bfgs") ,
  stmv_distance_statsgrid = 3, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** ),
  stmv_distance_prediction_max = 3 * 1.25 , # upper limit in distnace to predict upon (just over the grid size of statsgrid) .. in timeseries can become very slow so try to be small
  stmv_distance_scale = c( 25, 35, 45 ), #likely must be over 30km, so 50 +/- 20km, should likely match the setting in ~ line 256
  stmv_clusters = list( scale=rep("localhost", ncpus), interpolate=rep("localhost", ncpus) )
) #End passing of parameters


if (0) {

  p$stmv_global_modelformula = formula( paste(
    ' snowcrab.large.males_abundance',
    ' ~ s( t, k=3, bs="ts") + s( tsd, k=3, bs="ts") + s( degreedays, k=3, bs="ts") ',
    ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
    ' + s( substrate.grainsize, k=3, bs="ts")  ',
    ' + s(pca1, bs="ts") + s(pca2, bs="ts")  '
  ))

  p$stmv_global_modelformula = formula( paste(
    ' snowcrab.large.males_abundance',
    ' ~ s( t, k = 4, bs="ts") + s( tmax, k = 4, bs="ts") + s( degreedays, k = 4, bs="ts") ',
    ' + s( log(z), k=4, bs="ts") + s( log(dZ), k=4, bs="ts") + s( log(ddZ), k=4, bs="ts") ',
    ' + s( substrate.grainsize, k=4, bs="ts") + s(pca1, k=4, bs="ts") + s(pca2, k=4, bs="ts")   '
  ))

  p = stmv_variablelist(p=p)  # decompose into covariates, etc

  o = snowcrab_stmv(p=p, DS="input_data" )  # create fields for

  global_model = gam(
    formula=p$stmv_global_modelformula,
    family=p$stmv_global_family,
    data = o,
    weights=o$wt,
    optimizer= p$stmv_gam_optimizer,
    na.action="na.omit"
  )


}

#Run the following line if you want to use maptools rather than GADMTools for mapping coastline
# p$DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs", coastline_source="mapdata.coastPolygon" )'

# range( INP$snowcrab.large.males_abundance )
# [1]   14.3 6675.0

# o = snowcrab_stmv(p=p, DS="stmv_inputs" )  # create fields for

#--------------------------------------------------------------
# Run the process
#--------------------------------------------------------------

stmv( p=p, runmode=c("globalmodel", "interpolate" )  ) #  for a clean start

snowcrab_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids (if any)
snowcrab_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids (if any)
snowcrab_stmv( p=p, DS="complete.redo" )
snowcrab_stmv( p=p, DS="baseline.redo" )
snowcrab_stmv( p=p, DS="map.all" )

global_model = stmv_db( p=p, DS="global_model")
summary( global_model )

par(mar=c(1,1,1,1)) #change plot margins for Rstudio
plot(global_model)


# 2018 results
# Family: gaussian
# Link function: log
#
# Formula:
# snowcrab.large.males_abundance ~ s(t, k=3, bs="ts") + s(tsd,
#     k=3, bs="ts") + s(tmax, k=3, bs="ts") + s(degreedays,
#     k=3, bs="ts") + s(log(z), k=3, bs="ts") + s(log(dZ),
#     k=3, bs="ts") + s(log(ddZ), k=3, bs="ts") + s(substrate.grainsize,
#     k=3, bs="ts") + s(pca1, k=3, bs="ts") + s(pca2, k=3,
#     bs="ts")
#
# Parametric coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   6.7468     0.0204     330   <2e-16
#
# Approximate significance of smooth terms:
#                              edf Ref.df     F p-value
# s(t)                        1.91      2  73.0 < 2e-16
# s(tsd)                      2.00      2  14.0 8.3e-07
# s(tmax)                     1.91      2  11.7 1.4e-06
# s(degreedays)               1.36      2  45.4 < 2e-16
# s(log(z))                   1.91      2 174.7 < 2e-16
# s(log(dZ))                  1.93      2  20.9 3.8e-10
# s(log(ddZ))                 1.95      2  61.5 < 2e-16
# s(substrate.grainsize) 1.99      2  51.0 < 2e-16
# s(pca1)                     2.00      2 100.4 < 2e-16
# s(pca2)                     1.94      2  95.3 < 2e-16
#
# R-sq.(adj) =  0.232   Deviance explained = 23.4%
# GCV = 6198.2  Scale est. = 6182      n = 7640



# -------------------------------------------------------------------------------------
# STEP TWO commercial presence /absence
#--------------------------------------------------------------------------------------

# presence-absence
# this takes about 40 hrs ... and 5-6 GB /process
# year.assessment = 2018

p = bio.snowcrab::load.environment( year.assessment=year.assessment )
ncpus = parallel::detectCores()

p = snowcrab_stmv( p=p, DS="parameters",
  variables=list(Y="snowcrab.large.males_presence_absence"),
  selection=list(
    type = "presence_absence",
    biologicals = list(
      sex=0, # male
      mat=1, # do not use maturity status in groundfish data as it is suspect ..
      spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
      len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
      ranged_data="len"
    ),
    survey=list(
      data.source = c("snowcrab", "groundfish"),  # add groundfish data too
      drop.unreliable.zeros.groundfish.data=TRUE, # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
      yr = p$yrs      # time frame for comparison specified above
    )
  ),
  DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs" )',
  aegis_project_datasources = c("speciescomposition" ), # c("speciescomposition", "speciesarea", "sizespectrum", "condition", "metabolism", "biochem")
  stmv_global_family = binomial( link="logit" ),
  stmv_global_modelengine ="gam",
  stmv_local_modelengine = "twostep",
  stmv_twostep_time = "gam",
  stmv_twostep_space = "fft",
  stmv_fft_filter="matern",  #  matern, krige (very slow), lowpass, lowpass_matern
  stmv_gam_optimizer=c("outer", "bfgs") ,
  stmv_distance_statsgrid = 3, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** ),
  stmv_distance_prediction_max = 3 * 1.25 , # upper limit in distnace to predict upon (just over the grid size of statsgrid) .. in timeseries can become very slow so try to be small
  stmv_distance_scale = c( 25, 35, 45 ), #likely must be over 30km, so 50 +/- 20km, should likely match the setting in ~ line 256
  stmv_clusters = list( scale=rep("localhost", ncpus), interpolate=rep("localhost", ncpus) )
)


if (0) {

  p$stmv_global_modelformula = formula( paste(
    ' snowcrab.large.males_abundance',
    ' ~ s( t, k=3, bs="ts") + s( tsd, k=3, bs="ts") + s( degreedays, k=3, bs="ts") ',
    ' + s( t, tsd, degreedays, bs="ts")  ',
    ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
    ' + s( log(z), log(dZ), log(ddZ),  bs="ts") ',
    ' + s( substrate.grainsize, k=3, bs="ts")  ',
    ' + s(pca1, k=3, bs="ts") + s(pca2, k=3, bs="ts") + s(pca1, pca2, k=8, bs="ts")     '
  ))

  p$stmv_global_modelformula = formula( paste(
    ' snowcrab.large.males_abundance',
    ' ~ s( t, k = 4, bs="ts") + s( tmax, k = 4, bs="ts") + s( degreedays, k = 4, bs="ts") ',
    ' + s( log(z), k=4, bs="ts") + s( log(dZ), k=4, bs="ts") + s( log(ddZ), k=4, bs="ts") ',
    ' + s( substrate.grainsize, k=4, bs="ts") + s(pca1, k=4, bs="ts") + s(pca2, k=4, bs="ts")   '
  ))

  o = snowcrab_stmv(p=p, DS="input_data" )  # create fields for

  global_model = gam(
    formula=p$stmv_global_modelformula,
    family=p$stmv_global_family,
    data = o,
    weights=o$wt,
    optimizer= p$stmv_gam_optimizer,
    na.action="na.omit"
  )

}

# p$DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs", coastline_source="mapdata.coastPolygon" )'
# o = snowcrab_stmv(p=p, DS="stmv_inputs" )  # create fields for

stmv( p=p, runmode=c("globalmodel", "interpolate" ) ) # no global_model and force a clean restart

snowcrab_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids
snowcrab_stmv( p=p, DS="complete.redo" )
snowcrab_stmv( p=p, DS="baseline.redo" )
snowcrab_stmv( p=p, DS="map.all" )

global_model = stmv_db( p=p, DS="global_model")
summary( global_model )


par(mar=c(1,1,1,1)) #change plot margins for Rstudio
plot(global_model, all.terms=TRUE, trans=bio.snowcrab::inverse.logit, seWithMean=TRUE, jit=TRUE, rug=TRUE )


# Family: binomial
# Link function: logit
#
# Formula:
# snowcrab.large.males_presence_absence ~ s(t, k=3, bs="ts") +
#     s(tsd, k=3, bs="ts") + s(tmax, k=3, bs="ts") + s(degreedays,
#     k=3, bs="ts") + s(log(z), k=3, bs="ts") + s(log(dZ),
#     k=3, bs="ts") + s(log(ddZ), k=3, bs="ts") + s(substrate.grainsize,
#     k=3, bs="ts") + s(pca1, k=3, bs="ts") + s(pca2, k=3,
#     bs="ts")
#
# Parametric coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.0302     0.0349    0.87     0.39
#
# Approximate significance of smooth terms:
#                              edf Ref.df Chi.sq p-value
# s(t)                        1.99      2  336.7 < 2e-16
# s(tsd)                      2.00      2   35.3 2.0e-08
# s(tmax)                     1.98      2   59.9 7.0e-14
# s(degreedays)               1.99      2  188.4 < 2e-16
# s(log(z))                   2.00      2 1282.7 < 2e-16
# s(log(dZ))                  1.99      2   68.7 9.1e-16
# s(log(ddZ))                 1.98      2  224.4 < 2e-16
# s(substrate.grainsize) 1.98      2  394.4 < 2e-16
# s(pca1)                     2.00      2  932.6 < 2e-16
# s(pca2)                     2.00      2 1291.8 < 2e-16
#
# R-sq.(adj) =   0.62   Deviance explained = 56.2%
# UBRE = -0.58807  Scale est. = 1         n = 36031
#

# collect all predictions into a single file and return:
# year.assessment=2017

#------------------------------------------------------------------------------
# STEP THREE- Predict biomass by weighting +/- with probabilities
#------------------------------------------------------------------------------

p = bio.snowcrab::load.environment( year.assessment=year.assessment )
p = snowcrab_stmv( p=p, DS="parameters",
  variables = list(Y="snowcrab.large.males_abundance"),
  selection=list(
    type = "biomass",
    biologicals=list(
      sex=0, # male
      mat=1, # do not use maturity status in groundfish data as it is suspect ..
      spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
      len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
      ranged_data="len"
    ),
    survey=list(
      data.source = c("snowcrab"),
      yr = p$yrs,      # time frame for comparison specified above
      drop.unreliable.zeros.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
    )
  )
)

interpolation.db( DS="fishable.biomass.redo", p=p  ) # combine habitat and abundance info and map
K = interpolation.db( DS="fishable.biomass.timeseries", p=p  )
if (0){
  str(K)
  table.view( K )
  plot( total ~ yr, K[K$region=="cfanorth", ], type="b")
  plot( total ~ yr, K[K$region=="cfasouth", ], type="b")
  plot( total ~ yr, K[K$region=="cfa4x", ], type="b")
}

interpolation.db( DS="fishable.biomass.map", p=p  )

figure.timeseries.snowcrab.habitat(p=p) # /bio.data/bio.snowcrab/assessments/2016/timeseries/interpolated/snowcrab.habitat.sa.png

figure.timeseries.snowcrab.habitat.temperatures(p=p) # /bio.data/bio.snowcrab/assessments/2016/timeseries/interpolated/mean.bottom.temp.snowcrab.habitat.png





### --------- prediction success:
set = snowcrab_stmv(p=p, DS="input_data" )

S = set[ , c("plon", "plat") ]

ii = array_map( "xy->1", S, gridparams=p$gridparams )
bs = bathymetry.db(p=p, DS="baseline")
bb = array_map( "xy->1", bs, gridparams=p$gridparams )
im = match(  ii, bb )
it = match( set$yr, p$yrs )

bm = interpolation.db( DS="fishable.biomass", p=p  )
spred = bm$m[cbind(im, it)]  # approximate match (ignoring seasonality)

summary ( lm(log(spred)~log(snowcrab.large.males_abundance), data=set, na.actio="na.omit" ) )
plot(log(spred)~log(snowcrab.large.males_abundance), data=set )
cor(log(spred),log(set$snowcrab.large.males_abundance), use="complete.obs")

# Call:
# lm(formula = log(spred) ~ log(snowcrab.large.males_abundance),
#     data = set, na.action = "na.omit")
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -2.3039 -0.5507  0.0589  0.6056  2.1379
#
# Coefficients:
#                                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)                          5.17806    0.03824   135.4   <2e-16
# log(snowcrab.large.males_abundance)  0.19035    0.00599    31.8   <2e-16
#
# Residual standard error: 0.795 on 5030 degrees of freedom
#   (2223 observations deleted due to missingness)
# Multiple R-squared:  0.167,	Adjusted R-squared:  0.167
# F-statistic: 1.01e+03 on 1 and 5030 DF,  p-value: <2e-16
#
# R> plot(log(spred)~log(snowcrab.large.males_abundance), data=set )
# R> cor(log(spred),log(set$snowcrab.large.males_abundance), use="complete.obs")
# [1]  0.4269

# determine presence absence(Y) and weighting(wt)
#      set$weekno = trunc(set$julian / 365 * 52) + 1
#      set$dyear = trunc(set$julian / 365 ) + 1




## END
