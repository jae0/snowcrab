
  # --------------------------------- multispecies analysis ---------------
  # TODO: complete; not yet completed

  ordinate = function(x, varname, threshold=0.1, transf="log") {

    k = 1e8         # a large constant number to make xtabs work

    x = x[is.finite (x[,varname]),]
    x$id = paste(x$trip, x$set,  sep="~")
    m = xtabs( as.integer(x$totno*k) ~ as.factor(id) + as.factor(spec), data=x ) / k

    i = unique(c(which(rowSums(m)/dim(m)[1] < threshold ), which( rowSums(m)==0 ) ))
    j = unique(c(which(colSums(m)/dim(m)[1] < threshold ), which( colSums(m)==0 ) ))

      if (length(i) > 0 ) m = m[ -i , ]
      if (length(j) > 0 ) m = m[ , -j ]


      if (transf == "log" ) m = log10(m + 1)
      ord = cca( m )
      sp.sc = scores(ord)$species
      si.sc = scores(ord)$sites

      scores = data.frame( id=rownames(si.sc), ca1=as.numeric(si.sc[,1]), ca2=as.numeric(si.sc[,2]) )
      scores$id = as.character(scores$id)

  }



