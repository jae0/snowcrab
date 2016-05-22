

inflection.pt = function (expvar, pred, ste) {        
  # to determine the inflection point (day when eggcolour/CC changes) and associated se

  aa = tapply(pred, expvar, mean)
  me= approx( x= aa, y=sort(unique(expvar)), xout= 0.5)   #interpolate
  plot(sort(unique(expvar)), aa)              # use unique since there are multiple days that are the same.

  bb=(tapply(ste,expvar, mean))
  se=approx( y= bb, x=sort(unique(expvar)), me$y)
  plot(sort(unique(expvar)), bb)
  return ( list(me,se))     # combine in a list to return
}
   

