
    offset.time = function (x, ny=1) { # take data from previous year-class and year (i.e., previous data from the same cohort)
      dx = dim(x)
      if (!is.null(dx) ) {
        n.ageclass = dim(x)[1]
        n.years = dim(x)[2]
        i.years = c(1:n.years)
        i.y.offset = i.years - ny
        i.y.offset = i.y.offset[ i.y.offset >= 1 ]
        i.ageclass = c(1:n.ageclass)
        i.ac.offset = i.ageclass - ny
        i.ac.offset = i.ac.offset[ i.ac.offset >= 1 ]
        y = x * NA # initialise
        y[ (ny+1):n.ageclass, (ny+1):n.years] = x[i.ac.offset, i.y.offset]
      } else {
        n.ageclass = length(x)
        i.ageclass = c(1:n.ageclass)
        i.ac.offset = i.ageclass - ny
        i.ac.offset = i.ac.offset[ i.ac.offset >= 1 ]
        y = x * NA # initialise
        y[ (ny+1):n.ageclass ] = x[i.ac.offset]
      }

      return(y)
    }


