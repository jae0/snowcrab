
  predict.maturity.males = function(cw.mm, chela.mm) {
    # to ID maturity of males using mm
    # usage: predict.maturiyt.males(cw.mm=.., chela.mm=..)
    test.male = -25.32404 * log(cw.mm) + 19.775707 * log(chela.mm) + 56.649941
    mat = test.male * NA # initialise a vector with NA's of the right length
    mat[ which( test.male <= 0 ) ] = immature
    mat[ which( test.male > 0 )  ] = mature

    ii =which(cw.mm<60)
    mat[ii] = immature
    return(mat)
  }


