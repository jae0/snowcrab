
# odds and odds ratios

unadjusted.odds = function( .model, lc, not.factor ) {
    out = NULL
    model.formula = effects.labels( .model )

    for ( kf in  model.formula$effects ) { 
      
      cat("\n")

      print (kf )

      single.factor.model = as.formula( paste(as.character(  eval(model.formula$dep)), "~",  as.character(eval(kf) )) )
      
      xtab.formula = as.formula( paste("~ ",  as.character(eval(kf) )) )

      indata = remove.na( lc, single.factor.model, not.factor )

      xtab = as.data.frame.table( xtabs(xtab.formula, data=indata) )
      xt = data.frame(Nlevel = t(t(as.vector(xtab$Freq))))
      xt$Level = paste(kf, xtab[,1], sep="")

      logistic.results = glm( single.factor.model, data=indata, family=binomial )   
      
      res = logistic.display.local( logistic.results ) 
      
      res2 = merge(xt, res, by="Level", all=T, sort=F)
      res2 = res2[ order( res2$OddsRatio, na.last=F ) , ]

      res2$Covariate = kf
      res2$DepVariable = model.formula$dep
 
      out = rbind( out, res2)

      print (Anova(logistic.results) )
      cat("\n")
      
      print( res2 )
      cat("\n\n\n ------------------- \n")

    }
    return(out)

}
     
