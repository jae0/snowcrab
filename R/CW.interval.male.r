
  CW.interval.male = function(instar){
    i = instar - 3
    cw0 = exp(1.917564 + 0.298914*i)
    cw1 = exp(1.917564 + 0.298914*(i+1))
    return( c(cw0, cw1) )
  }


