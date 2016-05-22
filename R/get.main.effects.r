
get.main.effects = function( lm.model, dependent=F ) {
  model.statement = formula( lm.model )
  all.variables = all.vars(model.statement)
  dependent.var = all.variables[1] 
  main.effects = setdiff( all.variables, dependent.var )
  
  out = NULL
  out = ifelse (dependent, dependent.var, main.effects) 

  return(out)
}


