


# examine egg colour using a binomial GLM

  require (effects) # for "effect" for adjusted means
  require(car) # for "Anova"

    workdir = file.path( project.datadirectory( "bio.snowcrab"), "output" )  # linux
  	p = bio.snowcrab::load.environment()


# recreate the data files, only if the data needs to be refreshed
#  eggdata( data.loc=workdir, redo.data=T )

# load the data
#  sc = eggdata( data.loc=workdir, data.file="sc" )
  sdt = eggdata( data.loc=workdir, data.file="sdt" )


# choose the time variable resolution
  time.resolution = "weekly"
  sdt$timevalue = recode.time.maturity( sdt$julianday, time.resolution )
  yearoffset =  recode.time.maturity( 365, time.resolution )  # used when adding another year to the time value (see below)
  end.of.year = recode.time.maturity( 365*c(0,1,2,3), time.resolution )

  predictions.time = recode.time.maturity( c(1:730), time.resolution )
  predictions.eggpr =   c(0,1,2,3,4)
  predictions.temperature = round( quantile(sdt$t, probs=c(0.05, 0.5, 0.95), na.rm=T ), 2)


# ID various categories; p = primiparous; m = multiparous
  p12  = which( sdt$sex==2 & sdt$mat==1 & sdt$shell<=2 & sdt$eggcol %in% c(1,2) )   # n= 10,830
  p23  = which( sdt$sex==2 & sdt$mat==1 & sdt$shell<=2 & sdt$eggcol %in% c(2,3) )
  p34  = which( sdt$sex==2 & sdt$mat==1 & sdt$shell<=2 & sdt$eggcol %in% c(3,4) )
  m123 = which( sdt$sex==2 & sdt$mat==1 & sdt$shell>=3 & sdt$cw<100 )   # There were 3 outliers, cw>100 not likely!!
  m12  = which( sdt$sex==2 & sdt$mat==1 & sdt$shell>=3 & sdt$cw<100 & sdt$eggcol %in% c(1,2) )  #  n=32,803
  m23  = which( sdt$sex==2 & sdt$mat==1 & sdt$shell>=3 & sdt$cw<100 & sdt$eggcol %in% c(2,3) )  # n=19,898
  m34  = which( sdt$sex==2 & sdt$mat==1 & sdt$shell>=3 & sdt$cw<100 & sdt$eggcol %in% c(3,4) ) #  were only 2 females with eggcol=4




#################### Analysis finally begins ...

  # compare primi vs multi
  primi = sdt[ c(p12), ]
  primi$egg = 0
  primi$egg[ which(primi$eggcol==1)] = 1
  primi.binned = compute.proportions( primi, var="egg" )

 # compare primi vs multi
  multi = sdt[ c(m12), ]
  multi$egg = 0
  multi$egg[ which(multi$eggcol==1)] = 1
  multi.binned = compute.proportions( multi, var="egg" )


  # prim and mult 1-2:

   primi$stage = "P"
   multi$stage = "M"
   x=rbind( primi,multi)
  #Create a binary variable with shell condition, 0's for all other conditions, the basis for using the logit model
  x$egg = 0
  x$egg[ which(x$eggcol==1)] = 1

  model1 =  glm( egg ~   timevalue + stage + timevalue:stage , data=x, family=binomial()  )
  Anova(model1) # Get chi-sq and significance

  coef(model1)
  summary(model1)


  # ------------------


# Primiparous: Egg 1 to 2 transition
  x = sdt[ p12 , ]

  #Create a binary variable with shell condition, 0's for all other conditions, the basis for using the logit model
  x$egg = 0
  x$egg[ which(x$eggcol==1)] = 1

  model1 =  glm( egg ~   timevalue + t  + timevalue:t  , data=x, family=binomial()  )
  Anova(model1) # Get chi-sq and significance

  res = all.effects( model1, xlevels=list(t=predictions.temperature, timevalue=predictions.time), type="response" )
  plot(res, multiline=T, rescale.axis=F)

  coef(model1)
  summary(model1)


 # plot( res,  "timevalue:t", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.time))

#  plot( res,  "timevalue:eggPr", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.time))

#  plot( res,  "t:eggPr", type="response", multiline=T, rescale.axis=F, xlab="Egg fullness", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.eggpr) )


  to.plot = effects.data.export( res[["timevalue:t"]] )
  to.plot = effects.data.export( res[["timevalue:eggPr"]] )
  to.plot = effects.data.export( res[["t:eggPr"]] )

  binned.data = compute.proportions(x, var="egg")

  plot(to.plot$timevalue, to.plot$predictedPr, col="blue", pch=".")
  points(binned.data$timevalue, binned.data$fraction, col="red", pch=20)
  abline( v=end.of.year )


  export.to.clipboard( to.plot)
  export.to.clipboard( binned.data )

  # julian day at 50% transition of egg stage
  covariates = list( temp=predictions.temperature, eggpr=predictions.eggpr )  # associated with values of the other main effects where predictions are to be made
  tr = LDx( model1, p=0.5, error.method="2main.1interaction", covariates=covariates )



# ---------------------------
# Primiparous: Egg 2 to 3 transition
  x = sdt[ p23 , ]

  # ...  the following is not needed as there are no cases ... but kept to remind that it is conducted with multiparous females
  # move #3's btwn day 0-181 to the 'next year'
  # i = which( x$julianday <= 181 & x$julianday >= 0 & x$egg==0 )
  # x$timevalue[i] = x$timevalue[i] + yearoffset

  #Create a binary variable with shell condition, 0's for all other conditions, the basis for using the logit model
  x$egg = 0
  x$egg[ which(x$eggcol==3)] = 1

 i = which( x$julianday <= 181 & x$julianday >= 0 & x$eggcol ==3 )
 x$timevalue[i] =  x$timevalue[i] + yearoffset


  model1 =  glm( egg ~   timevalue +  t + timevalue:t, data=x, family=binomial()  )
  Anova(model1) # Get chi-sq and significance

  res = all.effects( model1, xlevels=list(t=predictions.temperature, timevalue=predictions.time), type="response" )

  plot(res, multiline=T, rescale.axis=F)

#  plot( res,  "timevalue", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.time) )

#  plot( res,  "eggPr:t", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.eggpr) )


  coef(model1)
  summary(model1)


  to.plot = effects.data.export( res[["timevalue:t"]] )
  to.plot = effects.data.export( res[["eggPr"]] )

  binned.data = compute.proportions(x, var="egg")

  plot(to.plot$timevalue, to.plot$predictedPr, col="blue", pch=".")
  points(binned.data$timevalue, binned.data$fraction, col="red", pch=20)
  abline( v=end.of.year )

  export.to.clipboard( to.plot)
  export.to.clipboard( binned.data)


  # julian day at 50% transition of egg stage
  covariates = list( temp=predictions.temperature,  eggpr=predictions.eggpr )  # associated with values of the other main effects where predictions are to be made
  tr = LDx( model1, p=0.5, error.method="3main.1interaction", covariates=covariates )



# ----------------------
# Multiparous: Egg 1 to 2 transition
  x = sdt[ m12 , ]

  #Create a binary variable with shell condition, 0's for all other conditions, the basis for using the logit model
  x$egg = 0
  x$egg[ which(x$eggcol==1)] = 1

  model1 =  glm( egg ~   timevalue + t + timevalue:t , data=x, family=binomial()  )
  Anova(model1) # Get chi-sq and significance

  res = all.effects( model1, xlevels=list(t=predictions.temperature, timevalue=predictions.time, eggPr=predictions.eggpr), type="response" )
  plot(res, multiline=T, rescale.axis=F)


 #  plot( res,  "timevalue:t", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.time))

 #  plot( res,  "timevalue:eggPr", type="response", multiline=T, rescale.axis=F, xlab="Egg fullness", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.time))

 #  plot( res,  "t:eggPr", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.temperature) )


  coef(model1)
  summary(model1)

  to.plot = effects.data.export( res[["timevalue"]] )
  to.plot = effects.data.export( res[["t"]] )
  to.plot = effects.data.export( res[["eggPr"]] )


  binned.data = compute.proportions(x, var="egg")

  plot(to.plot$timevalue, to.plot$predictedPr, col="blue", pch=".")
  points(binned.data$timevalue, binned.data$fraction, col="red", pch=20)
  abline( v=end.of.year )

  export.to.clipboard( to.plot)
  export.to.clipboard( binned.data)



  # julian day at 50% transition of egg stage
  covariates = list( temp=predictions.temperature, eggpr=predictions.eggpr )  # associated with values of the other main effects where predictions are to be made
  tr = LDx( model1, p=0.5,  error.method="2main.1interaction", covariates=covariates )



# ------------------------
# Multiparous: Egg 2 to 3 transition
  x = sdt[ m23 , ]

  # Create a binary variable with shell condition, 0's for all other conditions, the basis for using the logit model
  x$egg = 0
  x$egg[ which(x$eggcol==3)] = 1

  # 212 end of July, no BL eggs in Aug.
  i = which( x$julianday <= 181 & x$julianday >= 0 & x$eggcol==3 )
  x$timevalue[i] =  x$timevalue[i] + yearoffset


  model1 =  glm( egg ~   timevalue +  t , data=x, family=binomial()  )
  Anova(model1) # Get chi-sq and significance

  res = all.effects( model1, xlevels=list(t=predictions.temperature, timevalue=predictions.time), type="response" )

  plot(res, multiline=T, rescale.axis=F)

  # plot( res,  "timevalue", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: dark orange vs black", main="", xlim=range(predictions.time) )

  # plot( res,  "eggPr", type="response", multiline=T, rescale.axis=F, xlab="Egg fullness", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.eggpr))

  # plot( res,  "t", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: light orange vs. dark orange", main="", xlim=range(predictions.temperature) )


  coef(model1)
  summary(model1)

  to.plot = effects.data.export( res[["timevalue"]] )
  to.plot = effects.data.export( res[["t"]] )
  to.plot = effects.data.export( res[["eggPr"]] )

  binned.data = compute.proportions(x, var="egg")

  plot(to.plot$timevalue, to.plot$predictedPr, col="blue", pch=".")
  points(binned.data$timevalue, binned.data$fraction, col="red", pch=20)
  abline( v=end.of.year )

  export.to.clipboard( to.plot)
  export.to.clipboard( binned.data)


  # julian day at 50% transition of egg stage
  covariates = list( temp=predictions.temperature, eggpr=predictions.eggpr )  # associated with values of the other main effects where predictions are to be made
  tr = LDx( model1, p=0.5,  error.method="2main.0interaction", covariates=covariates )






# ---------------------------------------------
# Multiparous: Egg 3 to 4 transition
  x = sdt[ m34 , ]

  #Create a binary variable with shell condition, 0's for all other conditions, the basis for using the logit model
  x$egg = 0
  x$egg[ which(x$eggcol==4)] = 1

 # ? needs more work here ... ?
  # i = which( x$julianday <= 270 & x$julianday >= 0  & x$egg==0 ??? )
  # x$timevalue[i] =  x$timevalue[i] + yearoffset



  model1 =  glm( egg ~    timevalue + t + eggPr + timevalue:t + t:eggPr , data=x, family=binomial()  )
  Anova(model1) # Get chi-sq and significance

  res = all.effects( model1, xlevels=list(t=predictions.temperature, timevalue=predictions.time, eggPr=predictions.eggpr), type="response" )
  plot(res, multiline=T, rescale.axis=F)


 #  plot( res,  "timevalue:t", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: black vs coccoon", main="", xlim=range(predictions.time))

  coef(model1)
  summary(model1)

  to.plot = effects.data.export( res[["timevalue"]] )
  to.plot = effects.data.export( res[["t"]] )
  to.plot = effects.data.export( res[["eggPr"]] )

  binned.data = compute.proportions(x, var="egg")

  plot(to.plot$timevalue, to.plot$predictedPr, col="blue", pch=".")
  points(binned.data$timevalue, binned.data$fraction, col="red", pch=20)
  abline( v=end.of.year )

  export.to.clipboard( to.plot)
  export.to.clipboard( binned.data)


  # julian day at 50% transition of egg stage
  covariates = list( temp=predictions.temperature, eggpr=predictions.eggpr )  # associated with values of the other main effects where predictions are to be made
  tr = LDx( model1, p=0.5,  error.method="3main.2interaction.12.23", covariates=covariates )






# -------------------------------------
# Analysis of female maturity schedules

  x = sdt
  x$maturity = NA
  x$maturity[ which(x$mat=="2") ]  = 0 # immature
  x$maturity[ which(x$mat=="1") ]  = 1 # mature
  x = x[ is.finite( x$maturity) & x$sex==2, ]
  x$maturity = as.factor(x$maturity)


  newlymoulted.mature = which( x$maturity==1 & x$shell %in% c(1) )
  immature = which( x$maturity == 0 )
  p = sort( unique( c( immature, newlymoulted.mature ) ) )
  x = x[ p,]

  time.res = "daily"
  predictions.time = recode.time.maturity( c(-100:365), time.res )
  x$timevalue = recode.time.maturity( x$julianday, time.res )
  i = which( x$julianday > 180  )
  x$timevalue[i] =  x$timevalue[i] - recode.time.maturity( 365, time.res )

  model1 =  glm( maturity ~ timevalue + t + timevalue:t, data=x, family=binomial()  )
  Anova(model1) # Get chi-sq and significance

  res = all.effects( model1, xlevels=list( t=predictions.temperature, timevalue=predictions.time ), type="response" )
  plot( res,  "timevalue:t", type="response", multiline=T, rescale.axis=F, xlab="Days", ylab="Proportion: immature vs. mature", main="", xlim=range(predictions.time) )
  coef(model1)
  summary(model1)

  to.plot = effects.data.export( res[["timevalue:t"]] )
  to.plot = effects.data.export( res[["t"]] )


  binned.data = compute.proportions(x, var="maturity")

  plot(to.plot$timevalue, to.plot$predictedPr, col="blue", pch=".")
  points(binned.data$timevalue, binned.data$fraction, col="red", pch=20)
  abline( v=end.of.year )

  export.to.clipboard( to.plot)
  export.to.clipboard( binned.data )


  # julian day at 50% transition of egg stage
  tr = NULL
  for (temps in predictions.temperature ) {
    tr = rbind( tr, LDx( model1, p=0.5, error.method="2main.1interaction", covariates=temps ) )  # ######## must choose correct method depending on model specification

  }
  tr







