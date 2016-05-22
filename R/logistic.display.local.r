
logistic.display.local = function (logistic.model, alpha = 0.05, decimal = 5)  {

    model <- logistic.model
    
    if (class(model)[1] != "glm" | class(model)[2] != "lm" | 
        model$family$family != "binomial") {
        stop("model not from logistic regression")
    }
   
    s0 = Anova(model)
    s1 <- summary(model)

    orci <- as.data.frame(s1$coefficients)

    colnames(orci) <- c("OddsRatio", paste("Lower", 100 - 100 * alpha, 
        "CI", sep = ""), paste("Upper", 100 - 100 * alpha, "CI", 
        sep = ""), "P-value")

    orci[, 3] <- exp(orci[, 1] + qnorm(1 - alpha/2) * orci[, 2])
    orci[, 2] <- exp(orci[, 1] - qnorm(1 - alpha/2) * orci[, 2])
    orci[, 1] <- exp(orci[, 1])
  
    a <- orci[rownames(orci) != "(Intercept)", ]
    a$LogLikelihood = as.numeric(logLik(model)) 
    a$N = length(model$y)
    a$AIC = s1$aic
    a$deviance = s1$deviance
    a$Chisq.LogRatio = s0$"LR Chisq"
    a$Chisq.df = s0$Df
    a$Chisq.Pvalue = s0$"Pr(>Chisq)"

    a = round(a, decimal)
    
    a$Level = rownames(a)
    
    return(a)
}


