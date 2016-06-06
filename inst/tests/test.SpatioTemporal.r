  RLibrary( "geostatsinla", "geostatsp", "mgcv", "nlme4" )

  p = bio.snowcrab::initialise.local.environment()

  p$regions = c("cfa4x", "cfanorth","cfasouth" )
  p$vars.to.model = c("R0.mass")

  p$auxilliary.data = c(
            "t", "tmean", "tmean.cl", "tamp", "wmin",
            "z", "substrate.mean", "dZ", "ddZ",
            "ca1", "ca2",
            "nss.rsquared", "nss.shannon",
            "smr", "Ea", "A", "qm", "mass",
            "Z", "Npred" )

  p$habitat.threshold.quantile  = 0.05

  set0 = habitat.model.db( DS="basedata", p=p, v="R0.mass" )

  set = set0[ which(set0$yr %in% c( 2008:2013)) , ]



  require(SpatioTemporal)

  # FORMATTING REQUIRED for SpatioTemporal
  obs = data.frame(
    obs = set$R0.mass,
    date = set$yr + set$julian/365,
    ID = 1:nrow(set)
  )

  covars =
  st.covars =

  st  = createSTdata( obs=obs, covars=covars, transform.obs="log", SpatioTemporal=st.covars )


  coordinates( set ) = ~ plon + plat
  proj4string( set ) = CRS("+proj=utm +zone=20 +datum=WGS84")  # CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")


  # single year time slice
  locs = as.matrix( coordinates(set) )
  nData = nrow( set )

  bubble(set, "R0.mass")
  spplot(set, "R0.mass")





