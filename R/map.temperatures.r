
  map.temperatures = function( p, outdir, plottimes, conversions, init.files ) {
    tx = get.temperatures(DS="OE.saved")
    variables="t"
    tx$sa = 1/tx$err  # this mapping routine expects a weighting variable ("sa")
    p$block = F  # override ... the data are already blocked
    p$maskres      = "-S15k"
    p$interpres    = "-nb"
    p$tension      = "-T0.4"  # 1= harmonic surface
    gmt.map.variables(tx, p, variables=variables, plottimes=plottimes, basedir=outdir, conversions=conversions, init.files=init.files)
  }


