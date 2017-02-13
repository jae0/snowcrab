

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
  # mat=1, # maturity status in groundfish data is suspect
  spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
  len= c( 95, 200 )/10, #  mm -> cm ; indicators.db in cm
  drop.groundfish.data=TRUE # from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable 
)
p$lbm_local_modelengine = "twostep"
# p$lbm_twostep_space = "spatial.process"
# p$lbm_twostep_space = "fft"
# p$lbm_twostep_space = "tps"
p$lbm_twostep_space = "krige"
# p$lbm_twostep_space = "tps"
p$lbm_gam_optimizer=c("outer", "bfgs") 

p$lbm_global_family = gaussian(link="log")
p$lbm_local_family = gaussian(link="log")  # after logit transform by global model, it becomes gaussian (logit scale)
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
Link function: log 

Formula:
snowcrab.large.males_abundance ~ s(t, k = 3, bs = "ts") + s(tmean.climatology, 
    k = 3, bs = "ts") + s(tsd.climatology, k = 3, bs = "ts") + 
    s(log(z), k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + 
    s(log(ddZ), k = 3, bs = "ts") + s(log.substrate.grainsize, 
    k = 3, bs = "ts") + s(ca1, k = 3, bs = "ts") + s(ca2, k = 3, 
    bs = "ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  7.41427    0.02779   266.8   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                                 edf Ref.df       F  p-value    
s(t)                       1.544e+00      2   1.423 0.152650    
s(tmean.climatology)       1.106e+00      2  18.339 2.42e-10 ***
s(tsd.climatology)         7.818e-08      2   0.000 0.407220    
s(log(z))                  2.000e+00      2 215.450  < 2e-16 ***
s(log(dZ))                 3.065e-05      2   0.000 0.782840    
s(log(ddZ))                1.877e+00      2   6.597 0.000898 ***
s(log.substrate.grainsize) 1.731e+00      2  26.089 2.23e-13 ***
s(ca1)                     2.270e-04      2   0.000 0.217706    
s(ca2)                     1.998e+00      2 154.687  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.0999   Deviance explained = 10.2%
GCV = 7.8194e+06  Scale est. = 7.8012e+06  n = 4829



# -------------------------------------------------
# presence-absence
p = bio.snowcrab::load.environment( year.assessment=current.year )
p$selection=list( 
  name = "snowcrab.large.males_presence_absence",
  type = "presence_absence",
  sex=0, # male
  # mat=1, # maturity status in groundfish data is suspect
  spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
  len= c( 95, 200 )/10, #  mm -> cm ; indicators.db in cm
  drop.groundfish.data=TRUE # from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable 
)
p$lbm_local_modelengine = "twostep"
# p$lbm_twostep_space = "spatial.process"
# p$lbm_twostep_space = "fft"
# p$lbm_twostep_space = "tps"
# p$lbm_twostep_space = "krige"
p$lbm_twostep_space = "tps"
p$lbm_gam_optimizer=c("outer", "bfgs") 


p$lbm_global_family = binomial()
p$lbm_local_family = gaussian()  # after logit transform by global model, it becomes gaussian (logit scale)

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
    k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts") + s(log.substrate.grainsize, 
    k = 3, bs = "ts") + s(ca1, k = 3, bs = "ts") + s(ca2, k = 3, 
    bs = "ts")

Parametric coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  2.63028    0.02003   131.3   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                             edf Ref.df  Chi.sq  p-value    
s(t)                       1.996      2  218.78  < 2e-16 ***
s(tmean.climatology)       1.987      2   12.38  0.00195 ** 
s(tsd.climatology)         1.999      2  582.77  < 2e-16 ***
s(log(z))                  1.988      2 3703.88  < 2e-16 ***
s(log(dZ))                 1.999      2   62.51 2.55e-14 ***
s(log(ddZ))                1.992      2   28.01 7.82e-07 ***
s(log.substrate.grainsize) 1.991      2  414.41  < 2e-16 ***
s(ca1)                     1.993      2  464.80  < 2e-16 ***
s(ca2)                     1.995      2  469.86  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.398   Deviance explained =   38%
UBRE = -0.51619  Scale est. = 1         n = 61091





    # collect all results into a single file and return:
    K = interpolation.db( DS="interpolation.simulation", p=p  )
    table.view( K )

    figure.timeseries.errorbars( Pmeta, outdir=outdir, fname=paste(vv, rr, sep=".") )

    ### --------- prediction success:

    set = snowcrab.db( DS="set.complete" )
    set = set[ set$yr %in% p$yrs ,]
    set$total.landings.scaled = scale( set$total.landings, center=T, scale=T )
    set = presence.absence( set, "R0.mass", p$habitat.threshold.quantile )  # determine presence absence(Y) and weighting(wt)
#      set$weekno = floor(set$julian / 365 * 52) + 1
#      set$dyear = floor(set$julian / 365 ) + 1


  # update data summaries of the above results
    biomass.summary.db("complete.redo", p=p) #Uses the model results to create a habitat area expanded survey index
    biomass.summary.survey.db("complete.redo", p=p)#Uses average surface area from the past 5 years if a habitat area expanded surface area is not possible





