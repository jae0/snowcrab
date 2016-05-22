

  no.grid.logbook = function(x) {

    x$landings = x$pro_rated_slip_wt_lbs * 0.45359237  # convert to kg
    x$effort = x$num_of_traps
    x$cpue = x$landings / x$effort
    x$yr = x$year

    return(x)
  }


