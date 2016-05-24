

require(arm)
require(car)
require(effects)
require(splines)


  data.file = file.path(  project.datadirectory("snowcrab"), "R", "det_georef.rdata" )

	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
 


  load( data.file )
  det = det[ which(det$yr >= 1998),]  # data of limited spatial extent prior to 1998

  x = prep.maturity.data( det, sex="male", quantiles=c(0,0.25,0.5,0.75,1) )
  rm( det)

  .m0 = as.formula( "mat ~ cw*(yr + reg + t + z )" )
  .m0a = as.formula( "mat ~ cw*(yr + reg + bs(t, df=2) + z )" )
  .m0b = as.formula( "mat ~ cw*(yr + reg + bs(t, df=3) + z )" )
  .m0c = as.formula( "mat ~ cw*(yr + reg + bs(t,df=4) + z )" )
  .m0d = as.formula( "mat ~ cw*(yr + reg + bs(t,df=5) + z )" )
  .m1 = as.formula( "mat ~ cw*(yr + reg + tq + zq )" )
  .m2 = as.formula( "mat ~ cw*(yr + reg + tq + bs(z, df=3) )" )  # bs = B-spline from "splines"
  .m2a = as.formula( "mat ~ cw*(yr + reg + tq + bs(z, df=2) )" )  # bs = B-spline from "splines"
  .m2b = as.formula( "mat ~ cw*(yr + reg + tq + bs(z, df=4) )" )  # bs = B-spline from "splines"
  .m2c = as.formula( "mat ~ cw*(yr + reg + tq + bs(z, df=5) )" )  # bs = B-spline from "splines"
  .m3 = as.formula( "mat ~ cw*(yr + reg + t + bs(z, df=3) )" )  # bs = B-spline from "splines"
  
  .m0.glm = glm( .m0, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m0a.glm = glm( .m0a, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m0b.glm = glm( .m0b, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m0c.glm = glm( .m0c, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m0d.glm = glm( .m0d, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m1.glm = glm( .m1, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m2.glm = glm( .m2, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m2a.glm = glm( .m2a, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m2b.glm = glm( .m2b, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m2c.glm = glm( .m2c, data=x, family=binomial(link="logit"), na.action ="na.omit" )   
  .m3.glm = glm( .m3, data=x, family=binomial(link="logit"), na.action ="na.omit" )   

  anova( .m0.glm, .m1.glm )  # categorical seems better (.m1)
  anova( .m0.glm, .m3.glm )  # using a B-spline seems better (.m3)
  anova( .m1.glm, .m3.glm )  # categorical is better (.m1) 
  anova( .m1.glm, .m2.glm )  # mixed categorical with B-spline for depth seems best (.m2)
  anova( .m2.glm, .m2a.glm, .m2b.glm, .m2c.glm )  # check if number of nodes makes a difference .. .m2c seems best
  anova( .m0.glm, .m0a.glm, .m0b.glm, .m0c.glm, .m0d.glm )  # check if number of nodes makes a difference .. ..m0c seems best

  AIC(  .m0.glm, .m0a.glm, .m0b.glm, .m0c.glm, .m0d.glm , 
    .m1.glm,  .m2.glm,  .m2a.glm, .m2b.glm, .m2c.glm , .m3.glm )  #AIC agrees with anova (above)

  # final try: using splines for both temp and depth:

  .m4 = as.formula( "mat ~ cw*(yr + reg + bs(t, df=4) + bs(z, df=5) ) +bs( julian, df=4) + yr:reg" )
  .m4.glm = glm( .m4, data=x, family=binomial(link="logit"), na.action ="na.omit" )   

  anova( .m2c.glm , .m0c.glm, .m4.glm )  # .m4 is the best:: spline model in temp and depth
  AIC( .m2c.glm , .m0c.glm, .m4.glm )  # AIC agrees

 # finalised models:
  .m = .m4
  .m.glm = .m4.glm

  rm( .m0.glm, .m0a.glm, .m0b.glm, .m0c.glm, .m0d.glm, .m1.glm, .m2.glm, .m2a.glm, .m2b.glm, .m2c.glm, .m3.glm, .m4.glm)
  gc()

  # extract information
  Anova (.m.glm)
  AIC(.m.glm)
  summary( .m.glm ) 
  
  # various effects
  
    res.region = effect ( "reg", .m.glm ) # from effects
    plot( res.region, type="response" )    
    
    res.depth = effect ( "bs(z, df=5)", .m.glm ) # from effects
    plot( res.depth, type="response" )    

    res.temp = effect ( "bs(t, df=4)", .m.glm ) # from effects
    plot( res.temp, type="response" )    

    res.year = effect ( "yr:reg", .m4.glm ) # from effects
    plot( res.year, type="response" )    


    --- bring in ER and add as a factor and/or do a multi-level analysis with ER or simply a corrleation analysis with ER and YR:REG


    # effect of region:
    res = all.effects ( .m.glm ) # from effects
    plot( res, type="response", multiline=F, rescale.axis=F)    
   
    #  or simply
    
