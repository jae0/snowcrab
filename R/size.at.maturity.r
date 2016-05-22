
size.at.maturity = function(r) {
  
  require (doBy)

	outputvector = NULL

     res = coef(summary(r))
#      res2 = confint(r)
#      fin = cbind(res[,c("Estimate","Std. Error")],res2)
#      colnames(fin) = c("Estimate","Std. Error", colnames(res2))
      cw50 = dose.LD50(r, lambda=c(1,NA)) # lambda is a vector of model coef, with NA for inverse prediction variable.
      names(cw50) = c("cw50", "cw50lower", "cw50upper")
      olist = c("Estimate","Std. Error")
      r1 = res[1,olist] ; names(r1) = c("a0", "a0.se" )
      r2 = res[2,olist] ; names(r2) = c("a1", "a1.se" )
      r3 = c(length(i), r$deviance, r$aic); names(r3)=c("n", "deviance", "aic")
      r4 = c(floor(r$converged)); names(r4)="converged"

      outputvector = c( r1, r2, r3, r4, cw50)
  
}



