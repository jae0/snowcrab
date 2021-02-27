
year.assessment = 2019

p = bio.snowcrab::load.environment( year.assessment=year.assessment )

# 11 hrs with these settings
ncpus = parallel::detectCores()
stmv_clusters = list(
  scale=rep("localhost", ncpus),
  interpolate=rep("localhost", ncpus)
)

p = snowcrab_stmv(
  p=p,
  project_class="stmv",
  stmv_variables=list(Y="snowcrab.large.males_abundance"),
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
  stmv_gam_optimizer = c("outer", "bfgs") ,
  stmv_distance_statsgrid = 3, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** ),
  stmv_distance_scale = c( 25, 35, 45 ), #likely must be over 30km, so 50 +/- 20km, should likely match the setting in ~ line 256
  stmv_clusters = stmv_clusters
)


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

global_model = stmv_global_model( p=p, DS="global_model")
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

p = bio.snowcrab::snowcrab_parameters(
  p=p,
  project_class="stmv",
  stmv_variables=list(Y="snowcrab.large.males_presence_absence"),
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

global_model = stmv_global_model( p=p, DS="global_model")
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
p = snowcrab_stmv(
  p=p,
  project_class="stmv",
  stmv_variables = list(Y="snowcrab.large.males_abundance"),
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
bs = bathymetry_db(p=p, DS="baseline")
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
#      set$weekno = floor(set$julian / 365 * 52) + 1
#      set$dyear = floor(set$julian / 365 ) + 1




## END

  # figures and tables related to fishery indices

#Pick whichever year reference below is correct (most often year.assessment...-1)
if (!exists("year.assessment")) {
    year.assessment=lubridate::year(Sys.Date())
    year.assessment=lubridate::year(Sys.Date()) - 1
  }
  p = bio.snowcrab::load.environment( year.assessment=year.assessment )


# ------------------------------------------
# Time-series: all interpolated data estimated from interpolated analysis
# BZ R0.mass is likely the only variable required

  figure.interpolated.results( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ), alt.zero.y=T )

# ------------------------------------------
# Time-series: immature male numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.immature.male( p, outdir=file.path( p$annual.results, "timeseries", "interpolated" ) )

 # ------------------------------------------
# Time-series: immature male skip-moulter numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.immature.skipmoulter.male( p, outdir=file.path( p$annual.results,  "timeseries", "interpolated" ) )

# ------------------------------------------
# Time-series: mature CC12 male numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.mature.CC12.male( p, outdir=file.path( p$annual.results,  "timeseries", "interpolated" ) )

# ------------------------------------------
# Time-series: mature CC34 male numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.mature.CC34.male( p, outdir=file.path( p$annual.results, "timeseries", "interpolated" ) )

# ------------------------------------------
# Time-series: mature CC5 male numerical abundance estimates from interpolated analysis (broken down by instar)
  figure.mature.CC5.male( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ) )


# ------------------------------------------
# Time-series: Exploitation rate
  figure.exploitation.rate( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ) )
  figure.exploitation.rate( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ), CFA4X.exclude.year.assessment=F )

# ------------------------------------------
# Time-series: Generic plots of all interpolated data
  # figure.timeseries.interpolated( p, outdir=file.path( p$annual.results, "timeseries",  "interpolated" ) )

# ------------------------------------------
# Time-series: Habitat variations (surface area of snow crab habitat)

# p$model.type = "gam.full" # choose method for habitat model :
# p$model.type = "gam.simple" # choose method for habitat model :  no longer used

p$habitat.threshold.quantile = 0.05 # quantile at which to consider zero-valued abundance
p$prediction_dyear = 9/12 # predict for ~ Sept 1
figure.timeseries.snowcrab.habitat( p=p)

# ------------------------------------------
# Time-series: Habitat variations (mean temperature of snow crab habitat)
figure.timeseries.snowcrab.habitat.temperatures( p=p)

# ------------------------------------------
# Fecundity estimated indirectly via total number of females of primi and muli and applying mean egg production from allometric estimate
  fecundity.indirect( p=p,  outdir=file.path( p$annual.results, "timeseries", "survey" ) )

# ------------------------------------------
# Fecundity estimated directly via interpolation and individual-based fecundity estimate
  figure.timeseries.fecundity( p=p, outdir=file.path( p$annual.results, "timeseries",  "interpolated"  ) )




# Tables for the SSR/RESDOC
# ---------------------------------------------
  # Tables of CC from trawl surveys > 95mmCW

  outvars = c("yr", "region", "vars", "total", "lbound", "ubound")
  yy = c(1999:p$year.assessment)
  rr = c("cfasouth", "cfanorth", "cfa4x")
  vv = c( "totno.male.com.CC1", "totno.male.com.CC2", "totno.male.com.CC3",
          "totno.male.com.CC4", "totno.male.com.CC5" )

  p$runindex = list(y=yy, v=vv )

  L = interpolation.db( DS="interpolation.simulation", p=p )

  L$vars = as.character(L$vars)

  for (y in yy) {
  for (r in rr) {
    i = which( L$yr==y & L$region==r )
#    isum = sum(L$total[i], na.rm=T)
    L[i, c("total", "lbound", "ubound")] =round(  L[i, c("total", "lbound", "ubound")], 4)
  }}

  scale.factor = 1e2
  L[, c("total", "lbound", "ubound")] = L[, c("total", "lbound", "ubound")] * scale.factor
  L$vars = factor(x=L$vars, levels=vv, ordered=T)

  Mn = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfanorth"),], drop.unused.levels = T ) / scale.factor
  dimnames(Mn)$vars = c(1:5)

  Ms = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfasouth"),], drop.unused.levels = T ) / scale.factor
  dimnames(Ms)$vars = c(1:5)

  Mx = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfa4x"),], drop.unused.levels = T ) / scale.factor
  dimnames(Mx)$vars = c(1:5)

 (Mnp = round(Mn/rowSums(Mn)*100,2))
 (Msp = round(Ms/rowSums(Ms)*100,2))
 (Mxp = round(Mx/rowSums(Mx)*100,2))

  latex( Mnp, file="", title="", label="table.CC.north.trawl", rowlabel="Year", cgroup="Carapace condition", na.blank=T, caption="Carapace condition of crab larger than 95 mm CW (percent by number) over time for N-ENS from trawl surveys. The transition from a spring to a fall survey occurred in 2002/2003.")

  latex(Msp, file="", title="", label="table.CC.south.trawl", rowlabel="Year", cgroup="Carapace condition", na.blank=T, caption="Carapace condition of crab larger than 95 mm CW (percent by number) over time for S-ENS from trawl surveys. The transition from a spring to a fall survey occurred in 2002/2003.")

  latex(Mxp, file="", title="", label="table.CC.4x.trawl", rowlabel="Year", cgroup="Carapace condition", na.blank=T, caption="Carapace condition (percent by number) over time for CFA 4X from trawl surveys.")




# ---------------------------------------------
  # Tables of fishable biomass from trawl surveys


  outvars = c("yr", "region", "vars", "total", "lbound", "ubound")
  yy = c(1999:p$year.assessment)
  rr = c("cfasouth", "cfanorth", "cfa4x")
  vv = c( "totmass.male.com", "R0.mass", "R0.no", "R1.no", "R2.no", "R3.no", "totno.male.imm", "totno.male.com.CC1to2", "totno.male.com.CC3to4",
          "totno.male.com.CC5" )

  p$runindex = list(y=yy, v=vv )

  L = interpolation.db( DS="interpolation.simulation", p=p )
  # L = K[K$yr %in% yy & K$vars %in% vv & K$region %in% rr, outvars]
  L$vars = as.character(L$vars)

  for (y in yy) {
  for (r in rr) {
    i = which( L$yr==y & L$region==r )
#    isum = sum(L$total[i], na.rm=T)
    L[i, c("total", "lbound", "ubound")] =round(  L[i, c("total", "lbound", "ubound")], 4)
  }}


   scale.factor = 1
  L[, c("total", "lbound", "ubound")] = L[, c("total", "lbound", "ubound")] * scale.factor
  L$vars = factor(x=L$vars, levels=vv, ordered=T)

  Mn = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfanorth"),], drop.unused.levels = T ) / scale.factor
  dimnames(Mn)$vars = vv

  Ms = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfasouth"),], drop.unused.levels = T ) / scale.factor
  dimnames(Ms)$vars = vv

  Mx = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfa4x"),], drop.unused.levels = T ) / scale.factor
  dimnames(Mx)$vars = vv

  options(width=240)
  options(digits=1)
  Mn
  Ms
  Mx



  latex(Mn, file="", title="", label="table.stats.north.trawl", rowlabel="Year", cgroup="Various statistics", na.blank=T, caption="Various statistics of crab for N-ENS from trawl surveys. ")

  latex(Ms, file="", title="", label="table.stats.south.trawl", rowlabel="Year", cgroup="Various statistics", na.blank=T, caption="Various statistics for S-ENS from trawl surveys. ")

  latex(Mx, file="", title="", label="table.stats.4x.trawl", rowlabel="Year", cgroup="Various statistics", na.blank=T, caption="Various statistics for CFA 4X from trawl surveys.")
