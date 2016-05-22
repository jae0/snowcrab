
  fecundity.allometric = function( cw, method="st.marie" ) {
    
    if (method == "st.marie") {
      # this is from Saint-Marie(1993, unpublished) ... look for a more recent summary
      # est = cw ^2.725 * 10^0.7311  
    }
    if (method == "moncton.primiparous" ) {
      est = cw^2.94 * 10^(-0.7229)  # R^2 = 0.88
    }
    if (method == "moncton.multiparous" ) {
      est = cw^2.5811 * 10^0.0158   ## R^2 = 0.607
    }
    return( est )
  }

