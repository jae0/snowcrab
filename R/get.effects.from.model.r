get.effects.from.model = function( lm.model  ) {

    model.statement = formula( lm.model )
    all.variables = all.vars(model.statement)
    
    if (length( all.variables ) < 2) {
      print (" No results " )
      return ()
    }

    dependent.var = all.variables[1] 
    all.terms = attr( terms( model.statement), "term.labels" )
    

    # general stats (unadjusted medians/quantiles) 
    # for main effects that are significant
    
    # unadjusted effects 
    lm.data = lm.model$model  # internal copy of the data table
    results.unadjusted = list()
    for ( effect.variable in all.terms ) {
      print("---")
      print ( paste( effect.variable, " -- unadjusted: ") )
      results.unadjusted[[effect.variable]] =
extract.summary.data.for.one.variable( lm.data[,1],lm.data[,effect.variable])
} 
    lm.data$alldata = 1
    results.unadjusted[["alldata"]] =
extract.summary.data.for.one.variable( lm.data[,1], lm.data[,"alldata"]) 
results.adjusted = NULL
    for ( effect.variable in all.terms ) {
      print( paste( effect.variable, "..." ) )
      res = try( effect( effect.variable, lm.model ), silent=T )
      if (class(res)=="try-error") {
        print( "There was a problem extracting the effects from:" )
        print( effect.variable  )
        print( "And here is the error message:")
        print( res )
        print( "Here is the data and a summary of the data:")
        print( lm.model$model[[effect.variable]] )
        print( summary(lm.model$model[[effect.variable]]) )
        print( "try combining some levels? ..." )
        plot( lm.model$model[[effect.variable]]  ) 
      } else {
      med = as.vector(summary(res, typical="median")$effect)
      res.summary = summary(res)
      res.effect = res.summary$effect
      m = expand.grid(dimnames(res.effect))  
      if (!is.na(dim(res.effect)[2])) {
        n = paste(as.character(m[,1]), as.character(m[,2]), sep=".")
      } else {
        n = as.character(m[,1])
      }
      
     
      results.adjusted = rbind(results.adjusted, 
        cbind(dependent.var, effect.variable, n , exp(med)-1, 
          exp(as.vector(res.summary$effect))-1, 
          exp(as.vector(res.summary$lower))-1, 
          exp(as.vector(res.summary$upper))-1  ) )

filename = paste("adjusted.means", dependent.var, effect.variable,"png", sep=".") 
filename = gsub(":", "__", filename, fixed=T) 
      png( filename=filename )
       plot(effect(effect.variable, lm.model ))
      dev.off()

#      win.metafile(filename = graphfile, width = 7, height = 5)
#      plot(effect(effect.variable, lm.model ))
#      dev.off() 
     }
    }
    rownames(results.adjusted) = NULL
    colnames(results.adjusted) = c("interval", "effect", "level", "median", "mean", "95%CI.lower", "95%CI.upper")
    results.adjusted = as.data.frame(results.adjusted)
    for (i in 4:6) results.adjusted[,i] = as.numeric(as.character(results.adjusted[,i])) 
    print(" *.PNG Figures have been saved in your work directory" )
    
    return(list( results.unadjusted=results.unadjusted,
results.adjusted=results.adjusted)) 

}


