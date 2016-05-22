
  offset.ageclass = function (x, ny=1) { # take data from previous year-class and year (i.e., previous data from the same cohort)
    dx = dim(x)
    n.ageclass = dim(x)[1]
    i.ageclass = c(1:n.ageclass)
    i.ac.offset = i.ageclass - ny
    i.ac.offset = i.ac.offset[ i.ac.offset >= 1 ]
    y = x * NA # initialise
    y[ (ny+1):n.ageclass, , ] = x[i.ac.offset, , ]
    y [ !is.finite(y) ] = 0
      
    return(y)
  }



