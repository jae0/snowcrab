

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
# prep data for modelling and interpolation

snowcrab_lbm(p=p, DS="prediction.surface.redo" )  # create fields for snowcrab



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
p = bio.snowcrab::snowcrab.parameters( p=p, DS="lbm", varname=p$selection$name  )

# o = snowcrab_lbm(p=p, DS="lbm_inputs" )  # create fields for 
DATA='snowcrab_lbm( p=p, DS="lbm_inputs" )'

p = lbm( p=p, DATA=DATA, tasks=c("initiate", "globalmodel") ) # 30 min
#   p = lbm( p=p, tasks=c( "stage0" ) ) # serial mode
#   p = lbm( p=p, tasks=c( "continue" ) )    
p = lbm( p=p, tasks=c( "stage1" ) ) #  8 hrs 
p = lbm( p=p, tasks=c( "stage2" ) ) #   1 hrs
p = lbm( p=p, tasks=c( "save" ) )

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
snowcrab.large.males_abundance ~ s(yr) + s(dyear, k = 3, bs = "ts") + 
    s(yr, dyear, k = 36, bs = "ts") + s(ca1, bs = "ts") + s(t, 
    bs = "ts") + s(tmean, bs = "ts") + s(tamplitude, bs = "ts") + 
    s(log(z), bs = "ts") + s(log(dZ), bs = "ts") + s(log(ddZ), 
    bs = "ts") + s(log.substrate.grainsize, bs = "ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)   5.7329     0.1573   36.46   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                              edf Ref.df       F  p-value    
s(yr)                       1.000      1 108.941  < 2e-16 ***
s(dyear)                    2.000      2  93.129  < 2e-16 ***
s(yr,dyear)                32.988     33  14.904  < 2e-16 ***
s(ca1)                      8.766      9   5.175 2.89e-07 ***
s(t)                        8.400      9  14.439  < 2e-16 ***
s(tmean)                    8.510      9  13.077  < 2e-16 ***
s(tamplitude)               8.582      9   4.608 1.89e-06 ***
s(log(z))                   8.343      9  43.103  < 2e-16 ***
s(log(dZ))                  9.000      9  16.428  < 2e-16 ***
s(log(ddZ))                 8.665      9  28.501  < 2e-16 ***
s(log.substrate.grainsize)  8.139      9  17.750  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.239   Deviance explained = 24.8%
GCV = 9.264e+06  Scale est. = 9.1241e+06  n = 6977




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
p = bio.snowcrab::snowcrab.parameters( p=p, DS="lbm", varname=p$selection$name  )

# o = snowcrab_lbm(p=p, DS="lbm_inputs" )  # create fields for 
DATA='snowcrab_lbm( p=p, DS="lbm_inputs" )'

p = lbm( p=p, DATA=DATA, tasks=c("initiate", "globalmodel") ) # 30 min
#   p = lbm( p=p, tasks=c( "stage0" ) ) # serial mode
#   p = lbm( p=p, tasks=c( "continue" ) )    
p = lbm( p=p, tasks=c( "stage1" ) ) #  8 hrs 
p = lbm( p=p, tasks=c( "stage2" ) ) #   1 hrs
p = lbm( p=p, tasks=c( "save" ) )

p = make.list( list( yrs=p$yrs), Y=p )
parallel.run( snowcrab_lbm, p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_lbm( p=p, DS="lbm.stats.redo" ) # warp stats to other grids
snowcrab_lbm( p=p, DS="complete.redo" )
snowcrab_lbm( p=p, DS="baseline.redo" )
snowcrab_lbm( p=p, DS="map.all" )

global_model = lbm_db( p=p, DS="global_model") 
summary( global_model )
plot(global_model)








    # collect all results into a single file and return:
    K = interpolation.db( DS="interpolation.simulation", p=p  )
    table.view( K )

    figure.timeseries.errorbars( Pmeta, outdir=outdir, fname=paste(vv, rr, sep=".") )

    ### --------- prediction success:

    set = snowcrab.db( DS="set.complete" )
    set = set[ set$yr %in% p$years.to.model ,]
    set$total.landings.scaled = scale( set$total.landings, center=T, scale=T )
    set = presence.absence( set, "R0.mass", p$habitat.threshold.quantile )  # determine presence absence(Y) and weighting(wt)
#      set$weekno = floor(set$julian / 365 * 52) + 1
#      set$dyear = floor(set$julian / 365 ) + 1


  # update data summaries of the above results
    biomass.summary.db("complete.redo", p=p) #Uses the model results to create a habitat area expanded survey index
    biomass.summary.survey.db("complete.redo", p=p)#Uses average surface area from the past 5 years if a habitat area expanded surface area is not possible





