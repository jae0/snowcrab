
  # ----------------------------------------------------------------------------------------------
  # predict the maturity of crabs based upon determined data
  # ----------------------------------------------------------------------------------------------

  predictmaturity = function(x, method="DFA") {
    # codes
    # global parameters:: mature, immature, female, male is defined in load.snowcrab.environment 
    
    # basic biological limits
      x$mat[ which( x$cw<50 ) ]  = immature
      x$mat[ which( x$cw>150 ) ] = mature
    
    # males: order is important .. last is most imperative
      x$mat[ which( x$sex == male & x$shell %in% c(4,5)) ] = mature
      x$mat[ which( x$sex == male & x$cw < 50) ] = immature
      x$mat[ which( x$sex == male & x$cw > 150) ] = mature
      x$mat[ which( x$sex == male & x$chela < 5) ] = immature
      x$mat[ which( x$sex == male & x$chela > 33) ] = mature
      
    # female: sequence is important .. last is most imperative
      x$mat[ which( x$sex == female & x$shell %in% c(4,5)) ] = mature
      x$mat[ which( x$sex == female & x$cw < 30) ] = immature
      x$mat[ which( x$sex == female & x$cw > 80) ] = mature
      x$mat[ which( x$sex == female & x$abdomen < 20) ] = immature
      x$mat[ which( x$sex == female & x$abdomen > 50) ] = mature
      x$mat[ which( x$sex == female & x$gonad %in% c(1,2,3) )] = immature
      x$mat[ which( x$sex == female & x$eggPr == 0 )] = immature
      x$mat[ which( x$sex == female & x$eggPr %in% c(1:4) )] = mature # berried females
      x$mat[ which( x$sex == female & x$eggcol %in% c(1,2,3,4) )] = mature

    if (method=="DFA") {
      # males
      test.male = -25.32404 * log(x$cw) + 19.775707 * log(x$chela) + 56.649941
      x$mat[ which(x$sex==male & test.male <= 0 ) ] = immature
      x$mat[ which(x$sex==male & test.male > 0 ) ] = mature
    # females
      test.female = 16.422757 * log(x$cw) - 14.756163 * log(x$abdomen) - 14.898721
      x$mat[ which(x$sex==female  & test.female > 0 ) ] = immature
      x$mat[ which(x$sex==female  & test.female <= 0) ] = mature
    }
    
    if (method=="logistic.regression") {
      require(mgcv)
      
      x$rs = NA
      x$mat.predicted = NA

      WW = which( x$mat %in% c(immature, mature) )  
      MM = which( x$sex==male )
      FF = which( x$sex==female )
      SC = which( x$shell %in% 1:5 )
      cleandata_mm = which( x$cw > 0 & x$chela > 0 ) 
      cleandata_ff = which( x$cw > 0 & x$abdomen > 0 ) 

      x$mat = as.factor(x$mat)
      x$sex = as.factor(x$sex)
      #      x$shell = as.factor(x$shell)

      imm = intersect( intersect( intersect( MM, WW ), SC ), cleandata_mm )
      iff = intersect( intersect( intersect( FF, WW ), SC ), cleandata_ff )

      modm1 = glm( mat ~ log(cw) + log(chela) + shell, data=x[imm,], family=binomial(link="logit") )
      modf1 = glm( mat ~ log(cw) + log(abdomen) + shell, data=x[iff,], family=binomial(link="logit") )
      
      x$rs[imm]  = residuals( modm1 )
      x$rs[iff]  = residuals( modf1 )
      
      # reject too large a residual variability
      rs3 =  which( abs(x$rs) > 1.96 )  # reject too large a residual variability -- log-scaled ... residual =1 is large
      # these have potentially problematic data inputs in either cw and/or chela. etc
      # assume CW is OK and set all other data in the main dataset to uncertain/missing
      x$chela[rs3] = NA
      x$abdomen[rs3] = NA
      

      # update data vectors removing "bad" data and redo model
      cleandata_mm = which( x$cw > 0 & x$chela > 0 ) 
      cleandata_ff = which( x$cw > 0 & x$abdomen > 0 ) 
      imm = intersect( intersect( intersect( MM, WW ), SC ), cleandata_mm )
      iff = intersect( intersect( intersect( FF, WW ), SC ), cleandata_ff )

      modm1 = glm( mat ~ log(cw) + log(chela) + shell, data=x[imm,], family=binomial(link="logit") )
      modf1 = glm( mat ~ log(cw) + log(abdomen) + shell, data=x[iff,], family=binomial(link="logit") )
 

      # require( boot )  # inv.logit
      # plot( modm1, trans=inv.logit,  rug=T, jit=T, scale=0 )
      # plot( modf1, trans=inv.logit,  rug=T, jit=T, scale=0 )

      # identify inds without maturity indications
      WW = which( !(x$mat %in% c(immature, mature) )  )
      imm = intersect( intersect( intersect( MM, WW ), SC ), cleandata_mm )
      iff = intersect( intersect( intersect( FF, WW ), SC ), cleandata_ff )
 
      x$mat.predicted[imm] = ifelse( predict( modm1, x[imm,], type="response" ) < 0.5, immature, mature )
      x$mat.predicted[iff] = ifelse( predict( modf1, x[iff,], type="response" ) < 0.5, immature, mature )

      x$mat[imm] = x$mat.predicted[imm] 
      x$mat[iff] = x$mat.predicted[iff] 

      x$rs = NULL
      x$mat.predicted = NULL


    }
  
    return(x)
  }


