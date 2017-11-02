

  region.decompose = function(region) {
    area = region #

    # snowcrab surveys
    if (region=="cfa20") area = c("C20-I", "C20-O")
    if (region=="cfa21") area = c("C21")
    if (region=="cfa22") area = c("C22-I", "C22-O")
    if (region=="cfa23") area = c("C23A", "C23B", "C23C", "C23D", "C23SL")
    if (region=="cfa24") area = c("C24A", "C24B", "C24C", "C24D", "C24E", "C24SL")
    if (region=="cfa4x") area = "C24W"
    if (region=="cfanorth") area = c( "C20-I", "C20-O", "C21", "C22-I", "C22-O" )
    if (region=="cfasouth") area = c( "C23A", "C23B", "C23C", "C23D", "C23SL",
                                     "C24A", "C24B", "C24C", "C24D", "C24E", "C24SL" )
    if (region=="cfa24slope") area = "C24SL"
    if (region=="cfa23slope") area = "C23SL"
    if (region=="cfaslope") area = c("C24SL", "C23SL")
    if (region == "cfaall") area = c("C20-I", "C20-O", "C21", "C22-I", "C22-O",
                                     "C23A", "C23B", "C23C", "C23D", "C23SL",
                                     "C24A", "C24B", "C24C", "C24D", "C24E", "C24SL",
                                     "C24W" )
     # groundfish surveys
      Vn = as.character(c(440:442))   # these are alphanumeric codes
      Vs = as.character(c(443:452))   # these are alphanumeric codes
      W  = as.character(c(453:466))   # these are alphanumeric codes
      X  = as.character(c(470:495))   # these are alphanumeric codes
      Ge = c("5Z1","5Z2","5Z3","5Z4","5Z5","5Z6","5Z7","5Z8")

      if (region == "v") area = c(Vn, Vs)
      if (region == "w") area = W
      if (region == "x") area = X
      if (region == "ge") area = Ge
      if (region == "vw") area = c(Vn, Vs, W)
      if (region == "vwx") area = c(Vn, Vs, W, X)
      if (region == "all") area = c(Vn, Vs, W, X, Ge)

     
    return(area)
  }


