
  get.annual.timeseries = function( p, sex, outtype="yearclass", region="cfanorth") {
    
    # kriged data
    K = kriging.db( DS="UK.conditional.simulation.K.complete", p=p )
    K = factor2character(K, c("vars","region","vario.model"))

    K$vars = as.character(K$vars)
    K$region  = as.character(K$region)

    cl = make.classes (sex)
    yclass = cl$yclass
    cats = cl$cats
    varmap = cl$varmap
    yrs =  as.character( sort( unique( K$yr[K$yr>=1998] )) ) # years with reliable data

    ir = which( K$region %in% region )
    out = get.structured.data (yclass, yrs, cats, varmap, K[ir,], outvar="total")
    out.lbound = get.structured.data (yclass, yrs, cats, varmap, K[ir,], outvar="lbound")
    out.ubound = get.structured.data (yclass, yrs, cats, varmap, K[ir,], outvar="ubound")
    out.sd =  (out.ubound - out.lbound)  # this gives the 95% CI
    out.sd = out.sd / 4  # 1/2 range estimates 2 SD; 1/4 range est of 1 SD

    if (outtype=="yearclass") {
      OO = list(TS=out, TS.lb=out.lbound, TS.ub=out.ubound, TS.sd=out.sd)
    }
    if (outtype=="cohorts") {
      OO = get.cohort.data (yclass, yrs, cats, varmap, out, outvar="total")
    }

    return(OO)
  }


