

univariate.analysis = function( .model, lc, not.factor ) {
    out = NULL
    model.formula = effects.labels( .model )
    
    for ( kf in  model.formula$effects ) { 
      
      single.factor.model = as.formula( paste(as.character(  eval(model.formula$dep)), "~",  as.character(eval(kf) )) )

      lc.na.removed = remove.na(lc, single.factor.model, not.factor )

      results.lm = lm( single.factor.model, data=lc.na.removed )   
      
      res = as.data.frame( Anova( results.lm )[1,] )
      res$No = length(lc.na.removed [,1])
      res$Covariate = kf
      res$DepVariable = model.formula$dep
 
      out = rbind( out, res )

    }
    out = out[ order(out[,4]) , ]
    return(out)
}
            

