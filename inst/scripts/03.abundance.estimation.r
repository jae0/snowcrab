require(aegis.env)

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
#stmv_distance_prediction_fraction = 1, # stmv_distance_prediction = stmv_distance_statsgrid * XX ..this is a half window km (default is 0.75)


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
p = snowcrab_stmv( p=p, DS="parameters",
  variables=list(Y="snowcrab.large.males_abundance"),
  selection=list(
    type = "abundance",
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
  stmv_variogram_method = "gstat",
  stmv_distance_statsgrid = 3, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** ),
  stmv_distance_prediction_fraction = 1, # stmv_distance_prediction = stmv_distance_statsgrid * XX ..this is a half window km
  stmv_distance_scale = c( 25, 35, 45 ), #likely must be over 30km, so 50 +/- 20km, should likely match the setting in ~ line 256
  stmv_clusters = list( rep("localhost", 8), rep("localhost", 8), rep("localhost", 8) )  # no of cores used made explicit.. must be same length as "stmv_distance_scale"
) #End passing of parameters


if (0) {
  p$stmv_global_modelformula = formula( paste(
    ' snowcrab.large.males_abundance',
    ' ~ s( t, k = 3, bs = "ts") + s( tsd, k = 3, bs = "ts") + s( tmax, k = 3, bs = "ts") + s( degreedays, k = 3, bs = "ts") ',
    ' + s( tmean.climatology, k = 3, bs = "ts") + s( tsd.climatology, k = 3, bs = "ts")',
    ' + s( t.range, k=3, bs="ts") + s( b.range, k=3, bs="ts") ',
    ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
    ' + s( log(substrate.grainsize), k=3, bs="ts") + s(pca1, k=3, bs="ts") + s(pca2, k=3, bs="ts")   '
  ))
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
# snowcrab.large.males_abundance ~ s(t, k = 3, bs = "ts") + s(tsd,
#     k = 3, bs = "ts") + s(tmax, k = 3, bs = "ts") + s(degreedays,
#     k = 3, bs = "ts") + s(tmean.climatology, k = 3, bs = "ts") +
#     s(tsd.climatology, k = 3, bs = "ts") + s(t.range, k = 3,
#     bs = "ts") + s(b.range, k = 3, bs = "ts") + s(log(z), k = 3,
#     bs = "ts") + s(log(dZ), k = 3, bs = "ts") + s(log(ddZ), k = 3,
#     bs = "ts") + s(log(substrate.grainsize), k = 3, bs = "ts") +
#     s(pca1, k = 3, bs = "ts") + s(pca2, k = 3, bs = "ts")
#
# Parametric coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   6.6260     0.0244     271   <2e-16
#
# Approximate significance of smooth terms:
#                              edf Ref.df      F p-value
# s(t)                        1.90      2  60.26 < 2e-16
# s(tsd)                      1.80      2   3.40   0.023
# s(tmax)                     1.88      2  21.68 4.9e-11
# s(degreedays)               1.93      2  10.56 1.2e-05
# s(tmean.climatology)        1.85      2  33.19 2.6e-16
# s(tsd.climatology)          1.96      2  47.03 < 2e-16
# s(t.range)                  1.55      2   1.39   0.164
# s(b.range)                  2.00      2  13.77 9.7e-07
# s(log(z))                   1.70      2 158.80 < 2e-16
# s(log(dZ))                  1.95      2  31.10 1.4e-14
# s(log(ddZ))                 1.87      2  28.48 6.0e-14
# s(log(substrate.grainsize)) 1.99      2  38.53 < 2e-16
# s(pca1)                     1.99      2  91.15 < 2e-16
# s(pca2)                     1.98      2  81.12 < 2e-16
#
# R-sq.(adj) =  0.243   Deviance explained = 24.5%
# GCV = 6423.2  Scale est. = 6400.2    n = 7640



# -------------------------------------------------------------------------------------
# STEP TWO commercial presence /absence
#--------------------------------------------------------------------------------------

# presence-absence
# this takes about 40 hrs ... and 5-6 GB /process
# year.assessment = 2018

p = bio.snowcrab::load.environment( year.assessment=year.assessment )

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
  stmv_variogram_method = "gstat",
  stmv_distance_statsgrid = 3, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** ),
  stmv_distance_prediction_fraction = 1, # stmv_distance_prediction = stmv_distance_statsgrid * XX ..this is a half window km
  stmv_distance_scale = c( 25, 35, 45 ), #likely must be over 30km, so 50 +/- 20km, should likely match the setting in ~ line 256
  stmv_clusters = list( rep("localhost", 8), rep("localhost", 8), rep("localhost", 8) )  # no of cores used made explicit.. must be same length as )
)

if (0) {
  p$stmv_global_modelformula = formula( paste(
    ' snowcrab.large.males_presence_absence',
    ' ~ s( t, k = 3, bs = "ts") + s( tsd, k = 3, bs = "ts") + s( tmax, k = 3, bs = "ts") + s( degreedays, k = 3, bs = "ts") ',
    ' + s( tmean.climatology, k = 3, bs = "ts") + s( tsd.climatology, k = 3, bs = "ts")',
    ' + s( t.range, k=3, bs="ts") + s( b.range, k=3, bs="ts") ',
    ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
    ' + s( log(substrate.grainsize), k=3, bs="ts") + s(pca1, k=3, bs="ts") + s(pca2, k=3, bs="ts")   '
  ))
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
# snowcrab.large.males_presence_absence ~ s(t, k = 3, bs = "ts") +
#     s(tsd, k = 3, bs = "ts") + s(tmax, k = 3, bs = "ts") + s(degreedays,
#     k = 3, bs = "ts") + s(tmean.climatology, k = 3, bs = "ts") +
#     s(tsd.climatology, k = 3, bs = "ts") + s(log(z), k = 3, bs = "ts") +
#     s(log(dZ), k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts") +
#     s(log(substrate.grainsize), k = 3, bs = "ts") + s(pca1, k = 3,
#     bs = "ts") + s(pca2, k = 3, bs = "ts")
#
# Parametric coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -1.5210     0.0619   -24.6   <2e-16
#
# Approximate significance of smooth terms:
#                                  edf Ref.df Chi.sq p-value
# s(t)                        2.00e+00      2 1446.6 < 2e-16
# s(tsd)                      1.99e+00      2   79.2 < 2e-16
# s(tmax)                     1.98e+00      2   98.7 < 2e-16
# s(degreedays)               2.00e+00      2  324.7 < 2e-16
# s(tmean.climatology)        1.94e-07      2    0.0    0.51
# s(tsd.climatology)          1.21e+00      2   32.1 4.3e-09
# s(log(z))                   2.00e+00      2 2195.6 < 2e-16
# s(log(dZ))                  1.05e+00      2   19.4 3.3e-06
# s(log(ddZ))                 2.00e+00      2  273.0 < 2e-16
# s(log(substrate.grainsize)) 1.96e+00      2  211.8 < 2e-16
# s(pca1)                     2.00e+00      2  518.3 < 2e-16
# s(pca2)                     1.99e+00      2 1229.1 < 2e-16
#
# R-sq.(adj) =  0.622   Deviance explained = 55.8%
# UBRE = -0.56211  Scale est. = 1         n = 35868


# collect all predictions into a single file and return:
# year.assessment=2017

#------------------------------------------------------------------------------
# STEP THREE- Predict biomass by weighting +/- with probabilities
#------------------------------------------------------------------------------

p = bio.snowcrab::load.environment( year.assessment=year.assessment )
p = snowcrab_stmv( p=p, DS="parameters",
  variables = list(Y="snowcrab.large.males_abundance"),
  selection=list(
    type = "abundance",
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
interpolation.db( DS="fishable.biomass.map", p=p  )

K = interpolation.db( DS="fishable.biomass.timeseries", p=p  )


if (0){
  str(K)
  table.view( K )
  plot( total ~ yr, K[K$region=="cfanorth", ], type="b")
  plot( total ~ yr, K[K$region=="cfasouth", ], type="b")
  plot( total ~ yr, K[K$region=="cfa4x", ], type="b")
}

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
# [1] 0.409

# determine presence absence(Y) and weighting(wt)
#      set$weekno = floor(set$julian / 365 * 52) + 1
#      set$dyear = floor(set$julian / 365 ) + 1


### ______________ TESTING __________________

  # currently supported:
  # z = depth (m)
  # dZ = bottom slope (m/km)
  # ddZ = bottom curvature (m/km^2)
  # substrate.grainsize = mean grain size of bottom substrate (mm)
  # t = temperature (C) – subannual
  # tlb = temperature lower 95% bound (C) –subannual
  # tub = temperature upper 95% bound (C) –subannual
  # tmean = mean annual temperature
  # tsd = standard deviation of the mean annual temperature
  # tmin = minimum value of temperature in a given year – annual
  # tmax = maximum value of temperature in a given year – annual
  # tamplitude = amplitude of temperature swings in a year (tmax-tmin) – annual
  # degreedays = number of degree days in a given year – annual

p$variables$COV = c(
  "t", "tub", "tsd", "tmin", "tmax", "tamplitude", "degreedays",
  "tmean.climatology", "tsd.climatology",
  "z", "dZ", "ddZ",
  "substrate.grainsize", "pca1", "pca2"    )

u = snowcrab_stmv(p=p, DS="input_data", alldata=TRUE )
u$Y = u$snowcrab.large.males_abundance


qn = quantile( u$Y[u$Y > 0], probs=p$habitat.threshold.quantile, na.rm=TRUE )
u$Y[ u$Y < qn ] = qn/10

v = gam( Y ~ s( t, k = 3, bs = "ts") + s( tsd, k = 3, bs = "ts") + s( tmax, k = 3, bs = "ts") + s( degreedays, k = 3, bs = "ts")
    + s(log(z), k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts")
    + s(log(substrate.grainsize), k = 3, bs = "ts")
    + s(pca1, k = 3, bs = "ts") + s(pca2, k = 3, bs = "ts"),
      data = u,
      family=gaussian(link="log")
    )

summary(v)



u$Y = u$snowcrab.large.males_presence_absence

v = gam( Y ~ s( t, k = 3, bs = "ts") + s( tsd, k = 3, bs = "ts") + s( tmax, k = 3, bs = "ts") + s( degreedays, k = 3, bs = "ts")
    + s(log(z), k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts")
    + s(log(substrate.grainsize), k = 3, bs = "ts")
    + s(pca1, k = 3, bs = "ts") + s(pca2, k = 3, bs = "ts"),
      data = u,
      family=binomial(link="logit")
    )

summary(v)


## END
