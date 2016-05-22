format.crude.odds = function(r) {
    dec = 2
   r$Odd = paste( round(r$OddsRatio,dec), " (", round(r$Lower95CI,dec), " - ", round(r$Upper95CI,dec), " )", sep="")
   
   r$Level2 = r$Level
   for ( i in 1:nrow(res) ) {
     r$Level[i] = sub(r$Covariate[i], "",  r$Level2[i], fixed=T )
   }   
   
   output=r[, c("Covariate", "Level", "Nlevel", "Odd", "Chisq.Pvalue") ]
   return(output)
}


