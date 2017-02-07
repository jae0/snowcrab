

# --------------------------------------------------------------
#  Ensure the following scripts complete without error:
#  these external dependencies permit lookup of data, for snowcrab.db("complete.redo") :
#  system.file(package="bio.indicators", "scripts", "01.indicators.r")  
#  system.file(package="bio.indicators", "scripts", "02.*.r")  


# 1. Define some additional starting parameters for debugging
#    choose various over-rides: these are initially defined in parameters.r

p = bio.snowcrab::load.environment( year.assessment=2016 )

p$regions = c("cfa4x", "cfanorth","cfasouth" )

p$vars.to.model = c("R0.mass")
# p$vars.to.model = c("R0.mass",  "R1.no")
# p$vars.to.model = c("R0.mass", "R0.no", "R1.no", "totno.female.primiparous","totno.female.multiparous", "totno.female.berried", "fecundity","totno.female.imm", "totno.male.imm" )
# p$vars.to.model = c("R0.no", "R1.no", "totno.female.primiparous","totno.female.multiparous", "totno.female.berried", "fecundity","totno.female.imm", "totno.male.imm" )

 # p$years.to.model=2005:2012


# --------------------------------------------------------------
# using environmental data ... estimate/lookup missing environmental data .. (t,z)
logbook.db( DS ="fisheries.complete.redo", p=p )
snowcrab.db( DS ="set.complete.redo", p=p )


# -------------------------------------------------------------------------------------
# prep data for modelling and interpolation

snowcrab_lbm(p=p, DS="prediction.surface.redo" )  # create fields for snowcrab

selection=list( 
  name = "snowcrab.large.males_abundance",
  type = "abundance",
  sex=0, # male
  # mat=1, # maturity status in groundfish data is suspect
  spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
  len= c( bio.snowcrab::mb(8), 200)/10, #  mm -> cm ; indicators.db in cm
  drop.groundfish.data=TRUE # from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable 
)

p = bio.snowcrab::snowcrab.parameters( p=p, DS="lbm", varname=selection$name  )
snowcrab_lbm(p=p, DS="lbm_inputs", selection=selection )  # create fields for 

p = make.list( list( yrs=p$yrs), Y=p )

DATA='snowcrab_lbm( p=p, DS="lbm_inputs" )'

p = lbm( p=p, DATA=DATA, tasks=c("initiate", "globalmodel") ) # 5 min
#   p = lbm( p=p, tasks=c( "stage0" ) ) # serial mode
#   p = lbm( p=p, tasks=c( "continue" ) )    
p = lbm( p=p, tasks=c( "stage1" ) ) #  8 hrs 
p = lbm( p=p, tasks=c( "stage2" ) ) #   1 hrs
p = lbm( p=p, tasks=c( "save" ) )
parallel.run( snowcrab_lbm, p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_lbm( p=p, DS="lbm.stats.redo" ) # warp stats to other grids
snowcrab_lbm( p=p, DS="complete.redo" )
snowcrab_lbm( p=p, DS="baseline.redo" )
snowcrab_lbm( p=p, DS="map.all" )




selection=list( 
  name = "snowcrab.large.males_presence_absence",
  type = "presence_absence",
  sex=0, # male
  # mat=1, # maturity status in groundfish data is suspect
  spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
  len= c( bio.snowcrab::mb(8), 200)/10, #  mm -> cm ; indicators.db in cm
  drop.groundfish.data=TRUE # from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable 
)


p = bio.snowcrab::snowcrab.parameters( p=p, DS="lbm", varname=selection$name  )
snowcrab_lbm(p=p, DS="lbm_inputs", selection=selection )  # create fields for 

p = make.list( list( yrs=p$yrs), Y=p )

DATA='snowcrab_lbm( p=p, DS="lbm_inputs" )'

p = lbm( p=p, DATA=DATA, tasks=c("initiate", "globalmodel") ) # 5 min
#   p = lbm( p=p, tasks=c( "stage0" ) ) # serial mode
#   p = lbm( p=p, tasks=c( "continue" ) )    
p = lbm( p=p, tasks=c( "stage1" ) ) #  8 hrs 
p = lbm( p=p, tasks=c( "stage2" ) ) #   1 hrs
p = lbm( p=p, tasks=c( "save" ) )
parallel.run( snowcrab_lbm, p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_lbm( p=p, DS="lbm.stats.redo" ) # warp stats to other grids
snowcrab_lbm( p=p, DS="complete.redo" )
snowcrab_lbm( p=p, DS="baseline.redo" )
snowcrab_lbm( p=p, DS="map.all" )





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





