 
  # ---------------------------------------------------------------------
  # predict weights of individuals that do not have weights recorded from groundfish surveys
  # (for the moment, obtain weight-cw relationship from groundfish surveys and apply to cw from
  # snowcrab surveys
  # ------------------------------------------------------------------

  predictweights = function( Y, parameterisation="bio") {
    
    # fill in missing weights from historical data as they used regression to do this .. 
    # a slight variation is that the model chosen in 2007 includes shell condition as well

    if (parameterisation=="bio") {
      # develop regression relationships for each of: using cw, sex, maturity, shell, weighted by sa
      
      x = Y

      x$log.mass = log(x$mass)
      x$log.cw = log(x$cw)
      x$log.chela = log(x$chela)
      x$log.abdom =log(x$abdomen)
      x0 = x
      ii = which( x$sex%in% c(male, female) )
      x$sex = factor(x$sex, levels=c(male,female) )

      mod = glm( log.mass ~ log.cw + log.cw:sex, weights=sa, data=x, na.action='na.exclude')
      x$predictedmass = predict( mod, x )
      x$rstand =  abs( x$log.mass - x$predictedmass ) 
      x$rstand =   x$rstand  / sqrt(var(x$rstand, na.rm=T))

      jj = which(x$rstand  > 3 )
      x$log.mass [ jj ] = x$predictedmass [ jj ]
      kk = which( !is.finite( x$log.mass ) )
      x$log.mass [ kk ] = x$predictedmass [ kk ]

      # still a few NA's: fill in with a simpler model dependent only upon cw ~ 500 cases 
      mod = glm( log.mass ~ log.chela, weights=sa, data=x, na.action='na.exclude')
      x$predictedmass = predict( mod, x )
      kk = which( !is.finite( x$log.mass ) )
      x$log.mass [ kk ] = x$predictedmass [ kk ]
      
      # still a few NA's: fill in with a simpler model dependent only upon cw ~ 500 cases 
      mod = glm( log.mass ~ log.abdom, weights=sa, data=x, na.action='na.exclude' )
      x$predictedmass = predict( mod, x )
      kk = which( !is.finite( x$log.mass ) )
      x$log.mass [ kk ] = x$predictedmass [ kk ]
      
      Y$mass = exp(x$log.mass)

      return(Y)

# summary(mod);Anova(mod) from 2007 data

#Call:
#lm(formula = mass ~ cw + cw:sex + shell, subset = ii, weights = sa, 
#    na.action = "na.exclude")

#Residuals:
#      Min        1Q    Median        3Q       Max 
#-0.174287 -0.003357 -0.000181  0.003389  0.141542 

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -7.464833   0.024413 -305.78  < 2e-16 ***
#cw           2.890774   0.003796  761.55  < 2e-16 ***
#shell2       0.053130   0.017581    3.02   0.0025 ** 
#shell3       0.125508   0.017627    7.12  1.1e-12 ***
#shell4       0.154578   0.018308    8.44  < 2e-16 ***
#shell5       0.164819   0.039400    4.18  2.9e-05 ***
#cw:sex2     -0.023230   0.000649  -35.81  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#Residual standard error: 0.00838 on 13236 degrees of freedom
#  (36127 observations deleted due to missingness)
#Multiple R-Squared: 0.985,	Adjusted R-squared: 0.985 
#F-statistic: 1.46e+05 on 6 and 13236 DF,  p-value: <2e-16 

#Anova Table (Type II tests)

#Response: mass
#          Sum Sq    Df  F value Pr(>F)    
#cw         49.76     1 708761.9 <2e-16 ***
#shell       0.05     4    186.7 <2e-16 ***
#cw:sex      0.09     1   1282.6 <2e-16 ***
#Residuals   0.93 13236                    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 


#> summary(mod);Anova(mod)

#Call:
#lm(formula = mass ~ cw + cw:sex, data = newdata, weights = sa,     na.action = "na.exclude")

#Residuals:
#       Min         1Q     Median         3Q        Max 
#-0.0840392 -0.0021296 -0.0000281  0.0019451  0.1309121 

#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -7.7028530  0.0014793   -5207   <2e-16 ***
#cw           2.9630510  0.0003513    8435   <2e-16 ***
#cw:sex2     -0.0150139  0.0000779    -193   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#Residual standard error: 0.00379 on 205565 degrees of freedom
#  (1109 observations deleted due to missingness)
#Multiple R-Squared: 0.997,	Adjusted R-squared: 0.997 
#F-statistic: 3.96e+07 on 2 and 205565 DF,  p-value: <2e-16 

#Anova Table (Type II tests)

#Response: mass
#          Sum Sq     Df  F value Pr(>F)    
#cw        1135.4      1 79181818 <2e-16 ***
#cw:sex       0.5      1    37130 <2e-16 ***
#Residuals    2.9 205565                    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

      }
    
      if (parameterisation == "moncton") {
        m.soft = filter.class(x, "m.soft")
        m.hard = filter.class(x, "m.hard")
        x$mass[m.hard] = x$cw[m.hard]^ 3.098 *  0.0002665
        x$mass[m.soft] = x$cw[m.soft]^ 3.524 *  0.00002995

        #  use same equations as for male ... must fix this
        f.soft = filter.class(x, "f.soft")
        f.hard = filter.class(x, "f.hard")
        x$mass[f.hard] = x$cw[f.hard]^ 3.098 *  0.0002665
        x$mass[f.soft] = x$cw[f.soft]^ 3.524 *  0.00002995
      }

      return(x)
  }




