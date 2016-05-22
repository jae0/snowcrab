
  expand.weights = function (x, w){
    # expand a list of values and weights to allow brute force
    # calculation of means, frequency distributions, etc.
    weights = sort(unique(w))

    out = NULL
    for (i in weights) {
      a = x[ which(w == i) ]
      b = rep.int(a, i)
      out = c(out, b)
    }
    return(out)
  }

 
