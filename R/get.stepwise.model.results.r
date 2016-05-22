

get.stepwise.model.results = function( model.full ) {

    model.stepwise = stepAIC( model.full, direction="both" )
    
    print ("The stepwise selected model")
    model.stepwise$anova
    
    summary( model.stepwise )
    
    Anova ( model.stepwise )
    
    print ("The stepwise selected model coefficients")
    print(as.data.frame(coef( model.stepwise )))
    
    return( model.stepwise )
}


