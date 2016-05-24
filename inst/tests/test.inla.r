
# Comparison of spatial and geostatistical approaches of gstat, sp, goestatsinla and inla
#  ... and intial attempt at spatio-temporal approach as well with inla

  
  # NOTE ON Matern parameterization:
  #  nu is the smoothness parameter:
  #   if nu=0.5 : Cov(Dist) == exponential Cov/Cor function
  #   if nu -> Inf .. Cov(dist) -> Gaussian Cov/Cor Func
   
  # INLA: 
  #   sigma_x^2 = marginal variance
  #   sigma_e^2 = nugget variance .. measurement error
  #
  #   scale parameter ( rho ) = 1/phi_inla > 0
  #   range parameter phi_inla = 1/rho
  #   shape (or "smoothness") parameter (nu) > 0 .. order of the modified Bessel Function of second kind
  #   dist = euclidean distance (s_i - s_j) / 
  #   K_nu = modified Bessel function of the second kind of order nu 
  #
  #   Cor(dist) = 2^(1-nu) / Gamma( nu ) * ( rho * dist )^nu * K_nu( rho*dist )
  #   Sigma = Cov(dist) = sigma_x^2 * Cor(dist)
  #
  #   or if y is a latent field and x = observations .. i.e. e_i is obervation error or "nugget"
  #   y(s) = x(s) +  e_i ;; 
  #   # then the marginal covariance function of y is:
  #   Cov[y(s)] = sigma_e^2 I + Sigma
  #
  #   and "inla::range" = phi_inla*l ; l=length of a cell .. ie:
  #   RES$inla$summary.hyperpar["Range for space", "mode"] * xres(loaFit$raster) == RES$par$summary["range", "mode"])
  
  # geostatsp:
  #   Cor(dist) = 2^(1-nu) / Gamma(nu) * ( sqrt(8*nu) * dist / phi_gsp)^nu * K_nu(sqrt(8*nu) *dist/phi_gsp ) 
  #   range parameter : phi_gsp  ~= distance at which correlation is <0.14 and decaying regardless of the shape parameter nu
      mymatern = function(u, phi_gsp, nu) {
        # implementation of Matern used in geostatsp
        # u is distance
        # phi_gsp is 
        # besselK is second kind of order nu
        uscale = sqrt(8 * nu) * u/phi_gsp
        res = (1/(gamma(nu) * 2^(nu - 1))) * uscale^nu * besselK(uscale, nu)
        res[u == 0] = 1
        res
      }
  
  # Wikipedia and RandomFields::matern uses  (and calls nu "kappa" ):
  #   range parameter (phi_rf) = phi_gsp / 2 
  #   Cov(dist) = sigma^2 * (Gamma(nu) * 2^(nu-1) )^-1 * ( sqrt(2*nu) * dist/phi_rf ) ^nu * K_nu( sqrt(2*nu)*dist/phi_rf) 
 
  # geoR and RandomFields::whittle uses (and calls nu "kappa" ):
  #   cor = (2^(-(nu - 1))) / gamma(nu) * ((dist/phi_geoR)^nu) * besselK( x=dist/phi_geoR, nu=nu ) ;; besselK = third kind 
  #   ...  range parameter (phi_geoR) = phi_gsp / sqrt( 8*nu ) 
 
  # gstat (and calls nu "kappa" )
  #   phi_gstat is range parameter where 95% of total asymptotic variance is reached
  #   psill * ( 1.0 - (pow(2.0, -( nu - 1.0))/Gamma( nu )) * pow(dist / phi_gstat, nu ) * K_nu(dist / phi_gstat, nu, 1.0) )






##############
# ------------
# ------------
# ------------
##############

# USEFUL FUNCTIONS (for INLA)

  II = function(x) {x}  # dummy function to return itself .. used below

  # see documentation here: http://www.r-inla.org/examples/tutorials/tutorial-on-option-scale-model
  RLibrary( "INLA" )
  inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors



# DATA

  # convert snow crab data into a spatial data frame

  p= list()
	p$init.files = loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  p$init.files = c( p$init.files, loadfunctions( c( "spacetime", "utility", "parallel", "bathymetry", "temperature" ) ) )

  p$libs = RLibrary( "mgcv", "chron", "lattice", "lattice", "grid", "fields", "parallel", 
                         "sp", "INLA", "geostatsinla", "geostatsp", "raster"  ) 
  

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
  
  set0 = set0[ which(set0$yr %in% c( 2004:2013)) , ]
  
  set0$B = set0$R0.mass # geostatsinla does not like "." 's? 
  # set0$R0.mass = NULL
  
  fnca1 = ecdf( set0$ca1 ) 
  set0$qca1 = fnca1(set0$ca1) 
  set0$cca1 = as.numeric( cut( set0$qca1, breaks=c(0, 1/3, 2/3, 1), labels=c(1, 2, 3)  ))  # for use with categorical example
  
  locs0  = as.matrix( set0[,c("plon", "plat")] )
  
  nData0 = nrow( set0 )
  nyrs = length( unique( set0$yr) )
  set0$yrindex = set0$yr - min(set0$yr)  + 1  # nermeic encoding of year



# --------------------------
# --------------------------
# --------------------------

# simple kriging using sp/gstat 

  require(sp)

  set = set0[ set0$yr==2013, ]
  
  coordinates( set ) = ~ plon + plat 
  proj4string( set ) = CRS("+proj=utm +zone=20 +datum=WGS84")  # CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")


  # single year time slice
  locs = as.matrix( coordinates(set) )
  nData = nrow( set )

  bubble(set, "B")
  spplot(set, "B")

  bb = read.table( find.ecomod.gis( "snowcrab.boundingbox") )
  names(bb) =c("lon", "lat" )
  bb = lonlat2planar( bb, p$internal.projection )


  V = variogram( B ~ 1, data=set  )
  Vtot = var( set$B) 
  
  V$gamma = V$gamma / Vtot  # normalise it
  vmat  =  fit.variogram(V, vgm(psill=0.5, model="Mat", range=50, nugget=0.5, kappa=2 )) ## gstat's kappa is the Bessel function's "nu" smoothness parameter
  
  plot( V, model=vmat )

  Vrange = vmat$range[2] # 95% of total variance 
  Vpsill = vmat$psill[2] / (vmat$psill[1] + vmat$psill[2] ) 
  Zannual =  Vtot 

  Q <- gstat( id="simple", formula = B~1, data=set, model=vmat )
  
  set.grid = expand.grid( 
    plon=seq(from=p$corners$plon[1], to=p$corners$plon[2], length.out=200 ) ,
    plat=seq(from=p$corners$plat[1], to=p$corners$plat[2], length.out=200 ) 
  )
  coordinates( set.grid) = ~ plon + plat 
  proj4string( set.grid ) = CRS("+proj=utm +zone=20 +datum=WGS84")  # CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

  Qp <- predict( Q, set.grid )
  spplot( Qp["simple.pred"], main = "ordinary kriging predictions")
  
  mean( Qp@data$simple.pred )




# --------------------------
# --------------------------
# --------------------------


# equivalent Random Field representation using geostatsinla
 
  require( geostatsinla)

  fm =  formula( B ~ 1 ) 
  # rw2 is an ~ GAM
  # inla.group discretizes the data into smaller number of unique values for 
  # numerical efficiency and avoiding singular solutions
  
  ft = glgm( data=set, formula=fm, family="normal",
    cells=100, shape=2, buffer=100, 
    priorCI=list(sd = c(0.1, 10), range=c(10, 100)), 
    # control.compute=list(dic=TRUE, cpo=TRUE),
    # control.predictor=list(compute=T), # compute marginal distributions for each value of the linear predictor
    # control.inla=list(int.strategy = "grid", diff.logdens = 4, strategy = "laplace", npoints = 21), # better handle on tails
    #debug=TRUE, verbose=TRUE,    
    control.compute=list(dic=TRUE)
  )

  names(ft$parameter)
  plot(ft$parameter[["range"]]$posterior, type="l", col="green" )
  lines(ft$parameter[["range"]]$prior )

  plot(ft$parameter[["sd"]]$posterior, type="l", col="green" )
  lines(ft$parameter[["sd"]]$prior )

  names(ft$raster)
  plot(ft$raster[["predict.mean"]])
  plot(ft$raster[["random.mean"]])


  (ft$parameters$summary)
                 mean         sd 0.025quant  0.5quant 0.975quant      mode          kld
(Intercept)  1.000684  0.2788304  0.4260607  1.008548   1.530799  1.023375 2.958545e-11
range       53.841793 10.2431609 36.9611059 52.708591  77.058985 50.401980           NA
sdNugget     1.773049         NA  1.5883112  1.752210   1.942670        NA           NA
sd           1.385588         NA  1.0899911  1.382066   1.726755        NA           NA

  R = ft$inla
  (R$summary.fixed)  # intercept
  

  post.se = inla.tmarginal( function(x){ sqrt(1/x) }, R$marginals.hyperpar[[1]] ) # transforms precision to se scale
  inla.zmarginal( post.se ) # general stats on the marginal of hyper params  ... ie, Nugget .. 


  post.range.se = inla.tmarginal( function(x){ sqrt(1/x) }, R$marginals.hyperpar[["Range for space"]] ) 
  inla.zmarginal( post.range.se )

  fp = data.frame( rasterToPoints( ft$raster[["predict.mean"]] ) )
  names(fp) = c("plon", "plat", "S")

  region = "cfanorth"
  region = "cfasouth"

  i = filter.region.polygon(x=fp[,c("plon", "plat")], region=region, planar=T)
  fp = fp[i,]
  require(lattice)
  levelplot( B~plon+plat, data=fp, aspect="iso") 
  exc30 = excProb( ft, 1, nuggetInPrediction = TRUE) # conditional probabilities that  y >1     
    # nuggetInPrediction argument can be set to TRUE to compute probabilities of new
    # observations Y i exceeding a threashold, with FALSE specifying exceedance probabilities for
    # λ(s)

  require( geostatsp )

    # conditional simulation: [U |Y ] 
    oneSim = geostatsp::RFsimulate(
        model = ft$param, 
        data = ft$data["resid"], 
        err.model = ft$param["nugget"], 
        x = raster(extent( set),nrow = 100, ncol = 100)
    )


  # Prior 95% intervals for σ and φ are specified, 
  # glgm creates Gamma priors for the precision 1/σ^2 and a scaled range parameter φ/δ (with δ being the cell size) having the 95% intervals specified
  # inla requires priors to be continuous, but are otherwise unrestricted
  # also, priors are set for log precisions, with prior distributions available including the log-Gamma and Normal
  # Priors for the remaining parameters can be specified with inla arguments such as 
  #   control.fixed=list(prec.intercept=0.01) 


  # glgm returns: 
  #   inla (raw inla results)
  #   parameters - prior and posteriors; and "parameters$summary"
  #   raster (stack) - posterior means of random effects and fitted values on link scale g[lambda(s)] == raster[["random.mean"]]
  #     -- raster[["predict.invlogit"]] == posterior means of lamda(s)
 

  

# --------------------------
# --------------------------
# --------------------------

  # INLA directly using a Random Field representation .. step by step  .. space only ( take a single time slice)
      
    # mesh for computation
    
    # boundary domain
    M0.domain = inla.nonconvex.hull( locs0, convex=10, resolution=75 )
 
    M0 = inla.mesh.2d (
        loc=locs0, # locations of data points
        boundary=M0.domain, 
        offset=c( 10, 100 ),  # how much to extend in the c(inner, outer) domains
        max.edge=c( 10, 100 ),  # max size of a triange (in, out)
        min.angle=c(21),   # min angle (in, out)
        cutoff=10  # min distance allowed
    )       
    plot(M0, asp=1 ) # visualise mesh
    

    M = inla.mesh.2d (
        loc=locs, # locations of data points
        boundary=M0.domain, # locations of borders 
        offset=c( 10, 100 ),  # how much to extend in the c(inner, outer) domains
        max.edge=c( 20,100),  # max size of a triange (in, out)
        min.angle=c( 21 ),   # min angle (in, out)
        cutoff=5  # min distance allowed
    )                                   
    
    plot( M )
 
  # create a mesh onto which we can represent the Latent Random (Gaussian) Field -- all data years 
    

  ###
  ### ------------------------ COMMON PART now finished ------------
  ###



    
    
    # create and spde object with Matern paramaterizations
    S = inla.spde2.matern( mesh=M, alpha=2 ) # alpha is the matern/Bessel smoothness parameter (nu) default is 2 
 
    # simple model so this is not needed but to be consistent:
    i <- inla.spde.make.index('i', n.spde=S$n.spde )  
    # i = 1:S$n.spde  # alternatively this


    # create observation/prediction weight matrices for projection from latent field to observations
    A = inla.spde.make.A( mesh=M, loc=locs )


    # stack all data together to feed into inla .. format is sparse matrix-like 
    # effects encodes indices for i (the nodes of the GF) and m (the intercept terms)
    Z = inla.stack( 
      data = list( B=set@data$B),
      A = list( A, 1 ),
      effects = list( i=i, m=rep(1, nData) ),
      tag="estim"
    )

    # model fitting .. 0 = no intercept as it is encoded in the effects above
    R = inla( 
      B ~ 0 + m + f(i, model=S), 
      data=inla.stack.data(Z), 
      control.predictor=list( A=inla.stack.A(Z) )
    )

    names(R)

    (R$summary.fixed)   # intercept
    (R$summary.hyperpar) # dispersion parameters of the effects ( spatial i's, spatial m's )


    # posterior marginals of hyperparameters:
    # names(R$marginals.hyperpar)
    # [1] "Precision for the Gaussian observations"  .. nugget 
    # [2] "Precision for space"                    
    # [3] "Range for space"                        
    #
    # :
    post.se = inla.tmarginal( function(x){ sqrt(1/x) }, R$marginals.hyperpar[[1]] ) # transforms precision to se scale
    plot(post.se)
    inla.zmarginal( post.se ) # general stats on the marginal
    inla.emarginal( II, post.se )  # expected value of the marginal
    inla.qmarginal( c(0.025, 0.5, 0.975), post.se )  # quantiles of the marginal
    inla.hpdmarginal( 0.95, post.se )  # highest probability density ( 95% .. similar to the 95% credible (quantile) interval above)

    # posterior marginals of the latent field parameters: log(kappa) and log(tau)  --- With INLA: 1/kappa == range parameter 
    # tau is the variance parameter
    post.lfp = inla.spde2.result( R, "i", S, do.transf=TRUE ) 
    inla.emarginal( II, post.lfp$marginals.kappa[[1]] ) # 1/kappa gives "range"
    vv = inla.emarginal( II, post.lfp$marginals.variance.nominal[[1]] ) # (sigma_x) ^2
    sqrt(vv) ## sd_x  .. 
    
    inla.emarginal( II, post.lfp$marginals.range.nominal[[1]] )  # range at which corrlation falls to 0.14



    plot( R$marginals.fixed[["m"]], type="l", xlab="Intercept" ) 
    plot( R$summary.random[["i"]][,c("ID", "mean")], type="l", xlab="Coefficient of effect" )
    plot( R$marginals.hyper[[1]], type="l", xlab="Phi" )
    plot.default( inla.tmarginal( function(x) {1/exp(x)}, R$marginals.hyperpar[[  ]]), xlab="Kappa", type="l")



    # prediction of the random field: two ways .. with analysis or after analysis
    # (aside: mesh points are already predicted )

    newlocations = rbind(  c( 330, 4790), c(788, 5080) )
    Anew = inla.spde.make.A ( mesh=M, loc=newlocations )  # sparse matrix of the projection coefficients
    

    # A) simple projection, post-analysis  
    Anew %*% R$summary.random$i$mean

    # B) or using inla functionality    
    inla.mesh.project( inla.mesh.projector( M, loc=newlocations), R$summary.random$i$mean )

    # C) or using inla in a global analysis
     
      # simple model so this is not needed but to be consistent:
      i <- inla.spde.make.index('i', n.spde=S$n.spde )  
      # i = 1:S$n.spde  # alternatively this
    
    Znew = inla.stack( 
      data = list( B=NA),
      A = list( Anew ),
      effects = list( i=i, m=rep(1, nrow(newlocations) ) ),
      tag="spatialpredictions"
    )
    
    Z = inla.stack( Z, Znew )
   
    # refit the same model
    R = inla( 
      B ~ 0 + m + f(i, model=S), 
      data=inla.stack.data(Z), 
      control.predictor=list( A=inla.stack.A(Z), compute = TRUE ) ## new addition here to compute linear predictions but this is CPU expensive
    )

    ipredictions = inla.stack.index( Z, tag="spatialpredictions")
    R$summary.linear.predctor[ ipredictions$data ,]

   
    # Prediction onto a grid using projections (method A/B) is faster/computationaly more efficient .. method of choice
    pG = inla.mesh.projector( M, xlim=p$corners$plon, ylim=p$corners$plat, dims=c(p$nplons, p$nplats) )
    # pG = inla.mesh.projector( M, dims=c(p$nplons, p$nplats) )
    
    Pmean = inla.mesh.project( pG, R$summary.random$i$mean )  # posterior mean
    Psd = inla.mesh.project( pG, R$summary.random$i$s )       # posterior SD
    image(log(Pmean))

    
    # Prediction of the response variable (B) .. ie. imputation
    Zimp = inla.stack( 
      data=list(B=NA), 
      A=list( Anew, 1 ), # 1 is for the intercept
      effects=list( i=1:S$n.spde, m=rep(1, nrow(newlocations) ) ), # 1 is for the intercept
      tag="prediction.response"
    )

    Z = inla.stack( Z, Zimp )
    R = inla(  
      B ~ 0 + m + f(i, model=S), 
      data=inla.stack.data(Z), 
      control.predictor=list( A=inla.stack.A(Z), compute = TRUE ) ## new addition here to compute linear predictions but this is CPU expensive
    )

    iimputations = inla.stack.index( Z, tag="prediction.response")
    R$summary.fitted.values[ iimputations$data ,]
    oo = R$marginals.fitted.values[ iimputations$data ]
    inla.hpdmarginal( 0.95, oo[[2]] ) 
    inla.zmarginal( oo[[2]])


    # or more directly:
    R$summary.fix[1,1] + Anew %*% R$summary.random$i$mean  # for the intercept
    sqrt( 1^2 + R$summary.fix[1,2]^2 + drop( Anew %*% R$summary.random$i$sd^2 ) ) # for SE 

    ## re-run without computation of marginal distributions on the grid nor quantiles ..
    Zgrid = inla.stack( 
      data=list(B=NA), 
      A=list(pG$proj$A, 1),
      effects=list(i=1:S$n.spde, m=rep(1, prod(pG$lattice$dims)) ),
      tag="prediction.response.gridded"
    )

    Zall = inla.stack( Z, Zgrid) 
    
    Rall = inla( 
      B ~ 0 + m + f( i, model=S ), 
      data=inla.stack.data(Zall),
      control.predictor=list( A=inla.stack.A(Zall), compute=TRUE),
      quantiles=NULL,
      control.results=list(return.marginals.random=FALSE, return.marginals.predictor=FALSE )
    )


    igr = inla.stack.index( Zall, "prediction.response.gridded" )$data
    
    require( gridExtra)
    grid.arrange(
      levelplot( Pmean, col.regions=topo.colors(99), main="latent field means", xlab="", ylab="", scales=list(draw=FALSE)),
      levelplot( matrix( Rall$summary.fitted[igr,1], pG$lattice$dims[1]), xlab="", ylab="", main="response mean", col.regions=topo.colors(99), scales=list(draw=FALSE)),
    
      levelplot( Psd, col.regions=topo.colors(99), main="latent field SD", xlab="", ylab="", scales=list(draw=FALSE)),
      levelplot( matrix( Rall$summary.fitted[igr,2], pG$lattice$dims[1]), xlab="", ylab="", main="response SD", col.regions=topo.colors(99), scales=list(draw=FALSE)),
      nrows=2 
    )


    oo = inla.spde2.result( R, "i", S, do.transf=TRUE )
    pp = inla.tmarginal( function(x) 1/x, R$marginals.hyperpar [[1]] ) # sigma^2(e) .. nugget

    plot( R$marginals.fixed[[1]], type="l", xlab="Beta" )
    
    

    
  # ----------
  # --- assume B ~ Gamma distribution -- to make it positive definite
  # additional covariates .. ~ kriging with external drift
    
   
    # create observation/prediction weight matrices for projection from latent field to observations
    A = inla.spde.make.A( mesh=M, loc=locs )
    
    # create and spde object with Matern paramaterizations
    S = inla.spde2.matern( mesh=M, alpha=2 ) # alpha is the matern smoothness parameter .. 2 is the default
  

    # stack all data together to feed into inla .. format is sparse matrix-like 
    # effects encodes indices for i (the nodes of the GF) and m (the intercept terms)
    
    # simple model so this is not needed but to be consistent:
    i <- inla.spde.make.index('i', n.spde=S$n.spde )  
    # i = 1:S$n.spde  # alternatively this


    Z = inla.stack( 
      data = list( B=B ),
      A = list( A, 1 ),
      effects = list( 
        list(i=i ), # GRF
        data.frame( Intercept=1, gtmean=inla.group( set$tmean ), gca1=inla.group( set$ca1 ) )
      ) ,
      tag="data"
    )

    mf = B ~ 0 + Intercept + f(gtmean, model="rw1") + f(i, model=S) #dic= 760.1247 .. lower than without space .. better
    # mf = B ~ 0 + Intercept + f(gtmean, model="rw1")  #dic= 898.4635


    # model fitting .. 0 = no intercept as it is encoded in the effects above
    R = inla( 
      mf, 
      family="Gamma",
      # verbose=TRUE,
      data=inla.stack.data(Z), 
      control.predictor=list( A=inla.stack.A(Z), compute=TRUE ),
      control.compute=(list(dic=TRUE))
    )

    R$dic$dic

    # posterior of the intercept
    R$summary.fixed

    # precisions of the hyper-parameters
    R$summary.hyper[,c(1,2)]

    # variance and range of the spatial random effect .. marginal distributions
    r = inla.spde2.result( R, "i", S, do.transf=TRUE )  
    inla.emarginal( II, r$marginals.variance.nominal[[1]] ) # variance
    inla.emarginal( II, r$marginals.range.nominal[[1]] )    # (practical) range

    plot( R$marginals.fixed[["Intercept"]] , type="l" )
    plot( R$summary.random[["gtmean"]][, 1:2], type="l", ylab="Coef", xlab="tmean", ylim=range( R$summary.random[["gtmean"]][,c( 2,4,6) ] ) ) 
    lines(  R$summary.random[["gtmean"]][, c(1,4)], lty=2 ) # 95% CI
    lines(  R$summary.random[["gtmean"]][, c(1,6)], lty=2 )

  

    # create grid for output
    require (splancs) # load "inout"
    pG = inla.mesh.projector( M, xlim=p$corners$plon, ylim=p$corners$plat, dims=c(p$nplons, p$nplats) )
    inside = inout( pG$lattice$loc, M0.domain$loc ) 

    # random field (latent)
    xmean = inla.mesh.project( pG, R$summary.random$i$mean )
    xsd   = inla.mesh.project( pG, R$summary.random$i$sd )


 
  #---------------------
  # "Semi-continuous models" or delta model or heirarchical model 
  # separate model formulations for zero and non-zero values
    
    r_i == rainfall amount at location i

    z_i ~ Bernoulli(p_) == occurence variable (PA /absence .. 1 if r_>0 else 0 )
    y_i ~ Gamma(a_i, b_i) == amount of rain NA if r_i==0 else r_i

    logit( p_i ) = alpha_z + x_i :: alpha_z is intercept and x_i is a random effect modeled as a GRMF SPDE
    
    E[ y_i ] = mu_i = a_i/b_i   --- expected value
    V[ y_i ] = a_i / b_i^2 = mu_i^2 / phi    --- variance

    phi == precision parameter 

    log( mu_i) = alpha_y + beta_x * x_i

    alpha_y == intercept
    beta_x == scaling parameter to x_i, shared with the logit model 


    # now, ... using snow crab data:

    # AB of values greater than 0
    AB = set$B
    AB[ AB <= 0 ] = NA
    hist(AB, "fd") 
    mean(AB, na.rm=T)

    # PA/absence
    PA = set$Y
    hist(PA)
    mean(PA)


    # mesh M for computation
    plot(M, asp=1 ) # visualise mesh
    
    # matern SPDE representation object using M
    S = inla.spde2.matern( M, alpha=2 )
    
    # projection matrix A to translate from mesh nodes to data nodes
    A = inla.spde.make.A( mesh=M, loc=locs )

    # data stack for occurence (PA)
    Z.PA = inla.stack( 
      tag="PA",
      data=list( PA=PA, Y=cbind( PA, NA) ) ,
      A=list(A,1),
      effects=list( list(iPA=1:S$n.spde), list(b0_PA=rep(1,length(PA)))) 
    )


    # data stack for abudance
    Z.AB = inla.stack( 
      tag="AB",
      data=list( AB=AB, Y=cbind( NA, AB) ) ,  ### NOTE -- stacking in next col similar to the latent state-space modelling approach
      A=list(A,1),
      effects=list( list(i), list(b0_AB=rep(1,length(AB)))) 
    )

    # data stack with both 
    Z.all = inla.stack( Z.PA, Z.AB )


    # separate fits for each of the PA/absence and AB
    R.PA <- inla( 
      PA ~ 0 + b0_PA + f(iPA, model=S), 
      family='binomial', 
      data=inla.stack.data(Z.PA), 
      control.compute=list(dic=TRUE),
      control.predictor=list(A=inla.stack.A(Z.PA), compute=TRUE)
    )


    R.AB <- inla(
      AB ~ 0 + b0_AB + f(iAB, model=S), 
      family='gamma',  # log transf by default .. (?)
      data=inla.stack.data(Z.AB), 
      control.compute=list(dic=TRUE),
      control.predictor=list(A=inla.stack.A(Z.AB), compute=TRUE)
    )

    
    # now the joint model .. with the spatial component of PA having presidence
    R.PA <- inla( 
      Y ~ 0 + b0_PA + b0_AB + f( iPA, model=S) + f(iAB, copy='iPA', fixed=FALSE),
      family=c('binomial', 'gamma'),
      data=inla.stack.data(Z.all), control.compute=list(dic=TRUE),
      control.predictor=list(A=inla.stack.A(Z.all), compute=TRUE)
    )


    # and the joint model .. with the spatial component of AB having presidence
    R.AB_PA <- inla( 
      Y ~ 0 + b0_PA + b0_AB + f( iAB, model=S) + f(iPA, copy='iAB', fixed=FALSE),
      family=c('binomial', 'gamma'),
      data=inla.stack.data(Z.all), control.compute=list(dic=TRUE),
      control.predictor=list(A=inla.stack.A(Z.all), compute=TRUE)
    )

    # and the joint model .. with the spatial component of AB having presidence
    R.AB_PA.nospace <- inla( 
      Y ~ 0 + b0_PA + b0_AB ,
      family=c('binomial', 'gamma'),
      data=inla.stack.data(Z.all), control.compute=list(dic=TRUE),
      control.predictor=list(A=inla.stack.A(Z.all), compute=TRUE)
    )




    # posterior of the Beta_x (the scaling parameter for AB relative to the random spatial effect on PA/absence) 
    round( rbind( beta_x.pa=R.PA$summary.hy[4, ], beta_x.ap=R.AB_PA$summary.hy[4, ]), 4)

                   mean     sd 0.025quant 0.5quant 0.975quant   mode
        beta_x.pa 0.6183 0.0901     0.4432   0.6174     0.7972 0.6147
        beta_x.ap 1.4808 0.1757     1.1429   1.4777     1.8336 1.4683

        # so, yes .. spatial effects upon AB is correlated
        # .. i.e, if Beta_x is non-zero then there is significant influence of the PA-absence spatial field upon 
        #    spatial field of AB 

        # from Krainski & Lindgren 2014, p 52: .. values for β_x are different due that they are fitted taking into account two
        # factores. One is relatinship existence between the spatial effect on both linear predictor.
        # We infer this looking if β_x is less or greater than zero. The second factor is related to
        # the variance of the effect on the first equation. The value considering the two different 
        # orders are equal only if the variance of the effect on first equation is one. In other words,
        # the probability of rain occurrence is correlated with the amoutg of the rain, like the
        # prefencential sampling problem [Diggle et al., 2010


    # DIC comparisons (need to sum the correct ones as they are local to each observation):
    ip <- inla.stack.index( Z.all, tag='PA')$data
    dic.pa <- c( 
      dic.p = sum( R.PA$dic$local.dic[ip], na.rm=TRUE),
      dic.a = sum( R.PA$dic$local.dic[-ip], na.rm=TRUE)
    )
    
    dic.ap <- c( 
      dic.p = sum( R.AB_PA$dic$local.dic[ip], na.rm=TRUE),
      dic.a = sum( R.AB_PA$dic$local.dic[-ip], na.rm=TRUE)
    )

    # compare joint vs separate models
    rbind( sep=c( R.PA$dic$dic, R.AB$dic$dic), joint.pa=dic.pa, joint.ap=dic.ap)

                    dic.p    dic.a
        sep      421.0851 758.7886
        joint.pa 396.5557 717.9570
        joint.ap 399.7917 718.7502

        # joint models perform better .. pa seems "better" than ap 
  
    # compare influence of space
    out = data.frame( rbind( 
      dic.nospace = c(
        PA=sum( R.AB_PA.nospace$dic$local.dic[ip], na.rm=TRUE),  # PA/absence from the joint model 
        AB=sum( R.AB_PA.nospace$dic$local.dic[-ip], na.rm=TRUE)  # AB
      ), 
      dic.space=dic.ap
    ))
    out$total = rowSums( out) 
    out

                    PA AB    total
        dic.nospace 530.5609  897.9120 1428.473
        dic.space   399.7917  718.7502 1118.542
    # strong spatial effect upon model fits



    # posterior of the alpha_PA:: logit intercept term for the Bernoulli process 
   rownames( R.PA$summary.fix )
     round( rbind (
      sep = R.PA$summary.fix, 
      joint.pa=R.PA$summary.fix[1,],
      joint.ap=R.AB_PA$summary.fix[1,]
    ), 4)

                    mean     sd 0.025quant 0.5quant 0.975quant    mode kld
        b0_PA -0.4888 1.5108    -3.9812  -0.3799     2.5774 -0.2260   0
        joint.pa    -0.0385 0.6263    -1.4554   0.0168     1.0527  0.1184   0
        joint.ap    -0.0210 0.6001    -1.3806   0.0338     1.0189  0.1268   0

    # i.e., joint models have a more precise estimate of alpha_PA .. . but all are still ~ 0 valued
   names( R.PA$marginals.fix ) 
      sapply( list( 
        sep=R.PA$marginals.fix[[1]], 
        joint.pa=R.PA$marginals.fix[[1]],
        joint.ap=R.AB_PA$marginals.fix[[1]]
      ), 
      function(m) inla.pmarginal(0, m)  # cumulative distribution at zero
    )

              sep  joint.pa  joint.ap 
        0.6549014 0.4906642 0.4776905 


   # posterior of the alpha_AB:: logit intercept term for the Gamma process 
   rownames( R.AB$summary.fix )
    round( rbind (
      sep = R.AB$summary.fix, 
      joint.pa=R.PA$summary.fix[2,],
      joint.ap=R.AB_PA$summary.fix[2,]
    ), 4)

                        mean     sd 0.025quant 0.5quant 0.975quant    mode kld
        b0_AB -0.3612 0.5079    -1.5443  -0.3091     0.5471 -0.2109   0
        joint.pa     -0.8261 0.3844    -1.6797  -0.7965    -0.1385 -0.7452   0
        joint.ap     -0.8389 0.4029    -1.7415  -0.8080    -0.1186 -0.7556   0

        # similar estimates but joint models slightly more precise and means slightly more negative valued
        
        # and on their respective response scales
        c( PA=binomial(link='logit')$linkinv( R.PA$summary.fix[1,1]), 
           AB=exp( R.PA$summary.fix[2,1]) 
        )
        
         PA AB 
        0.4903791 0.4377397
        
        # vs observed values  ... SOME DIVERGENCE ...!
        c( PA=mean(PA, na.rm=TRUE), 
           AB=mean(AB, na.rm=TRUE)
        )

         PA AB 
        0.6495098 1.8587104 
        
  # posterior of the precision parameter of Gamma model
    rownames( R.PA$summary.hy )
    R.AB_PA$summary.hy[1, ]


    # Posterior marginal distributions of α_PA , α_AB , φ and β

    par(mfrow=c(2,2), mar=c(3,3,.5,.5), mgp=c(1.5,.5,0), las=1)
    plot( R.AB_PA$marginals.fix[[1]], type='l', ylab='Density', xlab=expression(alpha[PA]))
    plot( R.AB_PA$marginals.fix[[2]], type='l', ylab='Density', xlab=expression(alpha[AB]))
    plot( R.AB_PA$marginals.hy[[1]], type='l', ylab='Density', xlab=expression(phi))
    plot( R.AB_PA$marginals.hy[[4]], type='l', ylab='Density', xlab=expression(beta))


  # significance of the spatial effect x 
  # 1. quantiles: 
    ordx <- order( R.AB_PA$summary.random$iPA$mean)
    par(mar=c(3,3,0.5,0.5), mgp=c(1.5,0.5,0), las=1)
    plot( R.AB_PA$summary.random$iPA$mean[ordx], type='l', ylab='x', ylim=range( R.AB_PA$summary.random$iPA[, 4:6]))
    for (i in c(4,6)) lines( R.AB_PA$summary.random$iPA[ordx,i], lty=2)
    abline(h=0, lty=3)


  # parameters of the spatial random effect
  R.AB_PA.f <- inla.spde2.result( R.AB_PA, 'i', S, do.transf=TRUE)    
  
    # posterior of the the marginal variance
    inla.emarginal( function(x) x, R.AB_PA.f$marginals.variance.nominal[[1]] )
    plot( R.AB_PA.f$marginals.variance.nominal[[1]], type="l", xlab=expression(sigma[x]^2), ylab='Density' )

    # practical range ( at which corellation drops to 0.14)
    inla.emarginal( function(x) x, R.AB_PA.f$marginals.range.nominal[[1]] )
    plot( R.AB_PA.f$marginals.range.nominal[[1]], type="l", xlab='Practical range', ylab='Density' )

    # kappa
    plot.default( inla.tmarginal(function(x) exp(x), R.AB_PA$marginals.hy[[3]]), type='l', xlab=expression(kappa), ylab='Density')




  ### Space and time model 
  range(set0$yr)

  k = length( unique( set0$yr) )
   
  # create a mesh onto which we can represent the Latent Random (Gaussian) Field -- all data years 
  M0 = inla.mesh.2d (
        loc=locs0, # locations of data points
        boundary=M0.domain, 
        offset=c( 10, 100 ),  # how much to extend in the c(inner, outer) domains
        max.edge=c( 20,200),  # max size of a triange (in, out)
        min.angle=c(21),   # min angle (in, out)
        cutoff=5  # min distance allowed
    )       
  # create and spde object with Matern paramaterizations
  S = inla.spde2.matern( mesh=M0, alpha=2 ) # alpha is the matern/Bessel smoothness parameter (nu) default is 2 
    
  iset <- inla.spde.make.index('i', n.spde=S$n.spde, n.group=k)

  A <- inla.spde.make.A(mesh=M0, loc=locs0, group=set0$yr, n.group=k ) 

  Z = inla.stack( tag='stdata', data=list( Y=set0$B), A=list(A,1), effects=list(iset, w=set0$dyear ))

  # prior for temporal autoregression coefficient -- N( mean=0, variance=5) and initial value guess  
  h.spec <- list(theta=list(initial=0.7, param=c(0, 5)))

  # treating time as a categorical variable (factor) .. expand.factor.strategy = "inla" for control.fixed.argument
    
  formulae <- Y ~ 0 + w + f(i, model=S, group=i.group, control.group=list(model='ar1', hyper=h.spec) )

  res <- inla(formulae, data=inla.stack.data(Z), 
    control.predictor=list(compute=TRUE, A=inla.stack.A(Z)),
    control.family=list(initial=20, fixed=TRUE),
    control.inla=list(strategy='gaussian'), 
    control.fixed=list(expand.factor.strategy='inla') 
  )

  
  # summary of the covariate coefficients (together the observed mean for each covariate level)
  round(cbind(observed=tapply(set$B, set$dyear, mean), res$summary.fixed), 4)


  # summary for the posterior marginal distribution of the temporal correlation
  round(res$summary.hy[3,], 5)

  # marginals distributions for random field parameters on the user scale
  rf <- inla.spde2.result(res, 'i', attr(x.k, 'spde'), do.transf=TRUE)

    par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=2:0)

    plot(res$marginals.hyper[[3]], type='l', xlab=expression(beta), ylab='Density')
    abline(v=rho, col=2)
    
    plot(rf$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[x]), ylab='Density')
    abline(v=params[1], col=2) 

    plot(rf$marginals.kappa[[1]], type='l', xlab=expression(kappa), ylab='Density')
    abline(v=params[2], col=2)

    plot(rf$marginals.range.nominal[[1]], type='l', xlab='range nominal', ylab='Density')
    abline(v=sqrt(8)/params[2], col=2)


    
  # posterior random field

  idat <- inla.stack.index(effdat, 'stdata') 

  cor(set$B, res$summary.linear.predictor$mean[idat])


  # prediction for each time and visualize 
    stepsize <- 4*1/111
    nxy <- round(c(diff(range(coords[,1])), diff(range(coords[,2])))/stepsize)
    pG <- inla.mesh.projector(prmesh1, xlim=range(coords[,1]),
    ylim=range(coords[,2]), dims=nxy)

    xmean <- list()
    for (j in 1:k)
    xmean[[j]] <- inla.mesh.project(
    pG, res$summary.random$i$mean[iset$i.group==j])

    
    require(splancs)
    xy.in <- inout(pG$lattice$loc, cbind(PRborder[,1], PRborder[,2]))
    for (j in 1:k) xmean[[j]][!xy.in] <- NA

    require(gridExtra)
    do.call(function(...) grid.arrange(..., nrow=4),
    lapply(xmean, levelplot, xlab='', ylab='',
    col.regions=topo.colors(16), scale=list(draw=FALSE)))


  # validation
    #

vdat <- data.frame(y=as.vector(y), w=ccov,
time=rep(1:k, each=n),
xcoo=rep(coords[,1], k),
ycoo=rep(coords[,2], k))[-isel, ]
Aval <- inla.spde.make.A(prmesh1, loc=cbind(vdat$xcoo, vdat$ycoo),
group=vdat$time)
stval <- inla.stack(tag='stval', data=list(y=NA),
A=list(Aval,1), effects=list(iset, w=vdat$w))

# Now, we just use a full data stack to fit the model

stfull <- inla.stack(effdat, stval)
vres <- inla(formulae, data=inla.stack.data(stfull),
control.predictor=list(compute=TRUE, A=inla.stack.A(stfull)), 
control.family=list(initial=20, fixed=TRUE),
control.inla=list(strategy='gaussian'),
control.fixed=list(expand.factor.strategy='inla'))

# pred vs obs
ival <- inla.stack.index(stfull, 'stval')$data
plot(vres$summary.fitted.values$mean[ival], vdat$y, asp=1, xlab='Posterior mean', ylab='Observed')
abline(0:1, col=gray(.7))





  # ----


  tdb <- raster(ncols=100, nrows=100, crs=proj4string( set ) )
  tdb = rasterize( set )


  tdb = temperature.db( p=p, DS="temporal.interpolation.redo", yr=2000 ) 







##############
# ------------
# ------------
# ------------
##############


# Full spatio-temporal model with "Semicontinuous (Gamma, Bernoulli)" or delta-type response variable
    
# USEFUL FUNCTIONS (for INLA)


  II = function(x) {x}  # dummy function to return itself .. used below

  # see documentation here: http://www.r-inla.org/examples/tutorials/tutorial-on-option-scale-model
  RLibrary( "INLA" )
  inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors



# DATA

  # convert snow crab data into a spatial data frame

  p= list()
	p$init.files = loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  p$init.files = c( p$init.files, loadfunctions( c("spacetime", "utility", "parallel", "bathymetry", "temperature" ) ) )

  p$libs = RLibrary( "mgcv", "chron", "lattice", "lattice", "grid", "fields", "parallel", 
                         "sp", "INLA", "geostatsinla", "geostatsp", "raster"  ) 
  

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
  set0$B = set0$R0.mass # geostatsinla does not like "." 's? 
 

## DEBUG::
  set0 = set0[ which(set0$yr %in% c( 2010:2013)) , ]
  set0$plon = set0$plon
  set0$plat = set0$plat

  
  locs0  = as.matrix( set0[,c("plon", "plat")] )
  nData0 = nrow( set0 )
  nyrs = length( unique( set0$yr) )
  set0$yrindex = set0$yr - min(set0$yr)  + 1  # numeric encoding of year


  # boundary domain
  M0.domain = inla.nonconvex.hull( locs0, convex=10, resolution=75 )

  M0 = inla.mesh.2d (
      loc=locs0, # locations of data points
      boundary=M0.domain, 
      offset=c( 10, 100 ),  # how much to extend in the c(inner, outer) domains
      max.edge=c( 10, 100 ),  # max size of a triange (in, out)
      min.angle=c(21),   # min angle (in, out)
      cutoff=10 # min distance allowed
  )       
  plot(M0, asp=1 ) # visualise mesh
  
  # create a mesh onto which we can represent the Latent Random (Gaussian) Field -- all data years 
    


  # AB of values greater than 0
  effdat = data.frame( set0 )
  

  ## rename Y -variables as they cannot be the same in the effects list
  AB = effdat$B   
  AB[ AB <= 0 ] = NA

  # PA/absence
  PA = effdat$Y
  effdat$Y = NULL  # conflicts with "Y' below 

  # intercept estimates
  effdat$b0_PA=rep(1, nData0)
  effdat$b0_AB=rep(1, nData0)
 

  # SPDE components

  # matern representation using mesh M
    S0 = inla.spde2.matern( M0, alpha=2 ) # alpha=2 is exponential correlation function

    # indices of SPDE 
    # only time and mesh dependence ,, not location! 
    i <- inla.spde.make.index('i', n.spde=S0$n.spde, n.group=nyrs )  
    iPA <- inla.spde.make.index('iPA', n.spde=S0$n.spde, n.group=nyrs )  
    iAB <- inla.spde.make.index('iAB', n.spde=S0$n.spde, n.group=nyrs )  

    # projection matrix A to translate from mesh nodes to data nodes
    A = inla.spde.make.A( mesh=M0, loc=locs0, n.group=nyrs, group=effdat$yrindex )


  ## datastack == Y == ( PAAbsence, AB)

  # data stack for occurence (PA)
    Z.PA = inla.stack( 
      tag="presence_absence",
      data=list( PA=PA, Y=cbind( PA, NA)) ,
      A=list(A,1),
      effects=list( iPA=iPA, effdat ) 
    )


    # data stack for abundance
    Z.AB = inla.stack( 
      tag="abundance",
      data=list( AB=AB, Y=cbind( NA, AB)) ,  ### NOTE -- stacking in next col similar to the latent state-space modelling approach
      A=list(A,1),
      effects=list( iAB=iAB, effdat ) 
    )

    # data stack with both 
    Z.all = inla.stack( Z.PA, Z.AB )

    rm( PA, AB) 

  # to see specification of priors/hyperpriors ...    
    # ?inla.models

  # prior for temporal autoregression coefficient -- N( mean=0, variance=5) and initial value guess  

    # other options for inla:
    # ncpus = 22
    # ncpus = 8
    # num.threads=ncpus
    # control.compute=list(dic=TRUE),
    
    # Priors for the observations: use "control.family":
    #     control.family=list( hyper=list( prec=list( prior="loggamma", param=c(0.1, 0.1))))  # to assign a the hyperparameter (precision parameter) to a log Gamma prior with (shape) a=0.1, and (inverse-scale) b=0.1
    #     control.family=list( hyper=list( prec=list( prior="gaussian", param=c(0, 1))))  # to assign a the hyperparameter (precision parameter) to a gaussian prior with mean 0, and precision 1.
    theta.observations.presence = list( prec=list( prior="gaussian", param=c(0, 1/10 ))) # 100 .. sd=10
    theta.observations.abundance = list( prec=list( prior="gaussian", param=c(0, 1/10 )))
    theta.observations.presence_abundance = list( 
          prec.AB=list( prior="gaussian", param=c(0, 1/10)), 
          prec.PA=list( prior="gaussian", param=c(0, 1/10)) )


    # Priors for fixed effects (slopes or beta)
    # fixed effects are gaussian by default ~N( 0, 1/tau ) 
    # control.fixed = list(prec.intercept = 0.001, prec = 0.001)
    theta.beta.presence = list( mean = 0.001, prec = 1/10 )
    theta.beta.abundance = list( mean = 0.001, prec = 1/10 )
    theta.beta.presence_abundance = list( mean = 0, prec = 1/10 )


    # Priors
    # rw2 .. 
    #     d^2x_i ~ N( 0, tau^(-1) ) 
    #     theta = log( tau )
    #     i.e., no priors as E[d^2x_i]=0 .. only precision hyperparameter scale needs to be defined
    theta.rw2 = list( prec = list(prior="loggamma", param=c(1, 5e-5)) ) #default
    theta.rw2 = list( prec = list(prior="loggamma", param=c(1, 1e-2 )) ) 
    
    # ar1 .. 
    #     x_1 ~ N(0, (tau(1-rho^2))^(-1) ) 
    #     x_i = rho * x_(i-1) + epsilon_i .. epsilon_i ~ N(0,tau^(-1) )  ... for i 2... n
    #     abs( rho) < 1
    #     theta_1 = log(kappa)  ... where kappa == marginal precision = tau(1-rho^2) . ~logGamma()
    #     theta_2 = log( (1+rho) / (1-rho) ) ; ~N()
    #     prior theta is logit correlation = ( logit(cor)) , prec )
    theta.ar1 = list( rho=list(prior="normal", param=c( 0, 0.15) ) )  # defaults
    theta.ar1 = list( theta1=list(param=c(1, 0.1) ), rho=list( param=c( 0.9, 1/0.1 ) ) )  # N( mean=0.9, var=0.1) 
    

    # matern.2d ..
    #     Corr(d) = 2^(nu-1)*GammaFunc(nu))^(-1) * (kappa*d)^nu * BetaFunc_nu(kappa*d) ;;; alpha=nu+d/2
    #     r = (8*nu)^(1/2) * kappa^(-1) .. range
    #     nu = 1,2,3 ( the shape parmeter is defined for these integers for now ) 
    #     hyperparameters: theta = ( theta1=tau or 'log precison', theta2=r or 'log range')
    #     where 1/tau is the marginal variance "simga_x^2" 
    # defaults .. both on log scales
    theta.spde2 = list( theta1 = list( prior="loggamma", param=c(1, 5e-5 )),    # log( 1/ variance )
                        theta2 = list( prior="loggamma", param=c(1, 0.01 ) ) )  # log (range)
    theta.spde2 = list( theta = list(param=c( 1, 1/0.01 ) ) )  

    
    
   
    # ---------------------
    # PA only 
    fmla = ( PA ~ 0 + b0_PA 
      + f(iPA, model=S0, group=iPA.group, control.group=list(model='ar1', hyper=theta.ar1 ), diagonal=1e-2) 
      + f(tmean, model="rw2", hyper=theta.rw2, diagonal=1e-2 ) )  # spline with diagonal to prevent numerical singularity
    fmly = "binomial"
    
    Z = Z.PA
    
    R = inla( formula=fmla, family=fmly, 
             data=inla.stack.data(Z), 
             control.predictor=list(A=inla.stack.A(Z), compute=TRUE), 
             control.fixed=theta.beta.presence, 
             # control.inla = list(strategy="laplace", npoints=21 ),  # more points for tails (default is 9)
             # control.compute = list(cpo=TRUE, pit=TRUE ),  # cpo=conditional predictive ordinate .. leave one out measures of fit to id extreme values (p(y_i|y_{-i}) .. ie. posterior value; # PIT=probability Integral Transforms Pr( y_i {new} <= y_i | y_{-i} ) .. ie on pr scale
             verbose=TRUE )

    summary(R)
   
    # default model settings:
    # inla.model.properties(<name>, <section>)
    #
    # initial values are on internal scale !
    # priors are defined on internal scale

    
    plot(R, plot.fixed.effects=TRUE)  # Inla calls these "unstructured effect" .. default: Beta ~ N( 0, 0.0001 )  
    plot(R, plot.random.effects=TRUE)  # Inla calls these "structured effect" .. default: f() ~ N( 0, tau ); log(tau) ~ logGamma(1, 0.00005)
    plot(R, plot.hyperparameters=TRUE)  # precision scale
  
    graphics.off()
    
    # to report on user scale ... SD scale = tau^(-1/2)
    s <- inla.contrib.sd(R, nsamples=1000) # posterior "structured variability" -- all random effects combined? .. on SD scale
    (s$hyper) 
    hist(s$samples)
    plot(density(s$samples,bw=.1),xlab="sigma",main="")
 
    
    ---
 summary(R)
Call:
c("inla(formula = fmla, family = fmly, data = inla.stack.data(Z), ",  "    verbose = TRUE, control.predictor = list(A = inla.stack.A(Z), ",  "        compute = TRUE), control.fixed = theta.beta.presence)" )

Time used:
 Pre-processing    Running inla Post-processing           Total 
         0.4269       4663.3755          0.4144       4664.2168 

Fixed effects:
        mean     sd 0.025quant 0.5quant 0.975quant   mode kld
b0_PA 0.9157 0.1185     0.6833   0.9156     1.1483 0.9155   0

Random effects:
Name	  Model
 iPA   SPDE2 model 
tmean   RW2 model 

Model hyperparameters:
                    mean    sd      0.025quant 0.5quant 0.975quant mode   
Theta1 for iPA       7.4105  0.0005  7.4091     7.4105   7.4117     7.4105
Theta2 for iPA      -5.8854  0.0008 -5.8872    -5.8853  -5.8838    -5.8853
GroupRho for iPA     0.9904  0.0000  0.9904     0.9904   0.9904     0.9904
Precision for tmean  1.1818  0.0007  1.1801     1.1818   1.1833     1.1818

Expected number of effective parameters(std dev): 11.63(0.4246)
Number of equivalent replicates : 140.64 

Marginal Likelihood:  -14171.79 
Posterior marginal


    # --------------------------
    # AB only
    fmla = ( AB ~ -1 + b0_AB 
      + f(iAB, model=S0, group=iAB.group, control.group=list(model='ar1', hyper=theta.ar1 ), diagonal=1e-2 ) 
      + f(tmean, model="rw2", hyper=theta.rw2, diagonal=1e-2 ) ) 
    fmly = "gamma"
    
    Z = Z.AB
     
    R = inla( formula=fmla, family=fmly, 
             data=inla.stack.data(Z), 
             control.predictor=list(A=inla.stack.A(Z), compute=TRUE), 
             control.fixed=theta.beta.abundance, 
             # control.inla = list(strategy="laplace", npoints=21 ),  # more points for tails (default is 9)
             # control.compute = list(cpo=TRUE, pit=TRUE ),  # cpo=conditional predictive ordinate .. leave one out measures of fit to id extreme values (p(y_i|y_{-i}) .. ie. posterior value; # PIT=probability Integral Transforms Pr( y_i {new} <= y_i | y_{-i} ) .. ie on pr scale
             verbose=TRUE )

    summary(R)
   

Call:
c("inla(formula = fmla, family = fmly, data = inla.stack.data(Z), ",  "    verbose = TRUE, control.predictor = list(A = inla.stack.A(Z), ",  "        compute = TRUE), control.fixed = theta.beta.abundance)" )

Time used:
 Pre-processing    Running inla Post-processing           Total 
         0.4088      22879.1059          0.5884      22880.1030 

Fixed effects:
         mean     sd 0.025quant 0.5quant 0.975quant    mode kld
b0_AB -0.2402 0.1508    -0.5362  -0.2402     0.0556 -0.2402   0

Random effects:
Name	  Model
 iAB   SPDE2 model 
tmean   RW2 model 

Model hyperparameters:
                                               mean    sd      0.025quant 0.5quant 0.975quant mode   
Precision parameter for the Gamma observations  2.4456  0.0052  2.4342     2.4452   2.4568     2.4454
Theta1 for iAB                                  1.7071  0.0022  1.7011     1.7070   1.7128     1.7071
Theta2 for iAB                                 -3.1101  0.0026 -3.1170    -3.1104  -3.1035    -3.1100
GroupRho for iAB                                0.9464  0.0001  0.9461     0.9464   0.9467     0.9464
Precision for tmean                            18.4469  0.0566 18.2861    18.4398  18.6006    18.4478

Expected number of effective parameters(std dev): 311.46(0.354)
Number of equivalent replicates : 3.779 

Marginal Likelihood:  -14862.62 
Posterior marginals for linear predictor and fitted values computed



-- for gaussian prior for Y


summary(R)

Call:
c("inla(formula = fmla, family = fmly, data = inla.stack.data(Z), ",  "    verbose = TRUE, control.predictor = list(A = inla.stack.A(Z), ",  "        compute = TRUE), control.fixed = theta.beta.abundance)" )

Time used:
 Pre-processing    Running inla Post-processing           Total 
         0.5314        603.9707          0.8313        605.3335 

Fixed effects:
        mean     sd 0.025quant 0.5quant 0.975quant   mode kld
b0_AB 1.4335 0.1236     1.1875   1.4345     1.6739 1.4365   0

Random effects:
Name	  Model
 iAB   SPDE2 model 
tmean   RW2 model 

Model hyperparameters:
                                        mean     sd       0.025quant 0.5quant 0.975quant mode    
Precision for the Gaussian observations   0.5975   0.0462   0.5072     0.5978   0.6881     0.6005
Theta1 for iAB                            0.2985   0.1783  -0.0835     0.3128   0.6119     0.3578
Theta2 for iAB                           -1.9902   0.1630  -2.2763    -2.0033  -1.6416    -2.0441
GroupRho for iAB                          0.9260   0.0264   0.8676     0.9285   0.9692     0.9346
Precision for tmean                      98.0724 115.2660  10.7534    63.5590 393.9884    27.6759

Expected number of effective parameters(std dev): 316.23(38.14)
Number of equivalent replicates : 3.722 

Marginal Likelihood:  -16193.92 
Posterior marginals for linear predictor and fitted values computed

   
   
    # --------------------------
    # AB_PA 
    fmla = ( Y ~ 0 + b0_PA + b0_AB 
      + f( iAB, model=S0, group=iAB.group, control.group=list(model='ar1', hyper=theta.ar1 ), diagonal=1e-2 ) 
      + f( iPA, copy='iAB', fixed=FALSE, diagonal=1e-2) 
      + f( z, model="rw2", hyper=theta.rw2, diagonal=1e-2 ) 
      + f( tmean, model="rw2", hyper=theta.rw2, diagonal=1e-2 ) )
#    fmly = c("binomial", "gamma")
    fmly = c("binomial", "gaussian")
    Z = Z.all
    
    R = inla( formula=fmla, family=fmly, 
             data=inla.stack.data(Z), 
             control.predictor=list(A=inla.stack.A(Z), compute=FALSE ), 
             control.fixed=theta.beta.presence_abundance, 
             # control.family=theta.observations.presence_abundance,
             verbose=TRUE )
 


Call:
c("inla(formula = fmla, family = fmly, data = inla.stack.data(Z), ",  "    verbose = TRUE, control.predictor = list(A = inla.stack.A(Z), ",  "        compute = FALSE), control.fixed = theta.beta.presence_abundance)" )

Time used:
 Pre-processing    Running inla Post-processing           Total 
         0.4591      31990.3432          0.2343      31991.0366 

Fixed effects:
        mean     sd 0.025quant 0.5quant 0.975quant   mode kld
b0_PA -0.728 2.2376    -5.1211  -0.7281     3.6614 -0.728   0
b0_AB -0.728 2.2376    -5.1211  -0.7281     3.6614 -0.728   0

Random effects:
Name	  Model
 iAB   SPDE2 model 
tmean   RW2 model 
iPA   Copy 

Model hyperparameters:
                                                  mean    sd      0.025quant 0.5quant 0.975quant mode   
Precision parameter for the Gamma observations[2]  2.2441  0.0005  2.2427     2.2440   2.2460     2.2442
Theta1 for iAB                                     1.2642  0.0001  1.2638     1.2642   1.2650     1.2643
Theta2 for iAB                                    -3.1723  0.0002 -3.1729    -3.1723  -3.1715    -3.1723
GroupRho for iAB                                   0.7303  0.0000  0.7302     0.7303   0.7305     0.7303
Precision for tmean                               54.6048  0.1209 54.3004    54.5968  54.8510    54.6053
Beta for iPA                                       1.2010  0.0001  1.2006     1.2010   1.2017     1.2010

Expected number of effective parameters(std dev): 617.88(4.389)
Number of equivalent replicates : 4.553 

Marginal Likelihood:  -57392.30 


Call:
c("inla(formula = fmla, family = fmly, data = inla.stack.data(Z), ",  "    verbose = TRUE, control.predictor = list(A = inla.stack.A(Z), ",  "        compute = FALSE), control.fixed = theta.beta.presence_abundance)" )

Time used:
 Pre-processing    Running inla Post-processing           Total 
         0.4615      82722.7068          0.2350      82723.4032 

Fixed effects:
        mean     sd 0.025quant 0.5quant 0.975quant   mode kld
b0_PA 0.5322 2.2364    -3.8587   0.5321     4.9193 0.5322   0
b0_AB 0.5322 2.2364    -3.8587   0.5321     4.9193 0.5322   0

Random effects:
Name	  Model
 iAB   SPDE2 model 
z   RW2 model 
iPA   Copy 

Model hyperparameters:
                                                  mean    sd      0.025quant 0.5quant 0.975quant mode   
Precision parameter for the Gamma observations[2]  2.4084  0.0315  2.3400     2.4114   2.4633     2.4193
Theta1 for iAB                                     1.6740  0.0119  1.6527     1.6732   1.6992     1.6707
Theta2 for iAB                                    -3.2796  0.0025 -3.2828    -3.2798  -3.2734    -3.2791
GroupRho for iAB                                   0.9486  0.0001  0.9484     0.9486   0.9489     0.9486
Precision for z                                   11.8349  0.0145 11.8154    11.8337  11.8733    11.8381
Beta for iPA                                       0.3634  0.0167  0.3422     0.3591   0.4036     0.3492

Expected number of effective parameters(std dev): 329.68(5.025)
Number of equivalent replicates : 8.533 

Marginal Likelihood:  -49229.44 
   


Call:
c("inla(formula = fmla, family = fmly, data = inla.stack.data(Z), ",  "    verbose = TRUE, control.predictor = list(A = inla.stack.A(Z), ",  "        compute = FALSE), control.fixed = theta.beta.presence_abundance)" )

Time used:
 Pre-processing    Running inla Post-processing           Total 
         0.4999      63567.7499          0.3120      63568.5618 

Fixed effects:
        mean     sd 0.025quant 0.5quant 0.975quant   mode kld
b0_PA 0.5668 2.2376    -3.8264   0.5667     4.9563 0.5668   0
b0_AB 0.5668 2.2376    -3.8264   0.5667     4.9563 0.5668   0

Random effects:
Name	  Model
 iAB   SPDE2 model 
z   RW2 model 
tmean   RW2 model 
iPA   Copy 

Model hyperparameters:
                                           mean    sd      0.025quant 0.5quant 0.975quant mode   
Precision for the Gaussian observations[2]  3.3288  0.0000  3.3288     3.3288   3.3379     3.3289
Theta1 for iAB                             -0.0935  0.0000 -0.0936    -0.0935  -0.0909    -0.0935
Theta2 for iAB                             -4.3087  0.0000 -4.3087    -4.3087  -4.3061    -4.3087
GroupRho for iAB                            0.6442  0.0000  0.6442     0.6442   0.6450     0.6442
Precision for z                            60.4302  0.0006 60.4287    60.4300  60.5896    60.4304
Precision for tmean                        52.9654  0.0005 52.9641    52.9653  53.1052    52.9656
Beta for iPA                                0.9544  0.0000  0.9544     0.9544   0.9571     0.9544

Expected number of effective parameters(std dev): 1109.92(0.00)
Number of equivalent replicates : 2.534 

Marginal Likelihood:  -39584.48 
> 


    # --------------------------
    # NO space 
    fmla = ( Y ~ -1 + b0_PA + b0_AB 
      + f( yr, model="rw2", hyper=theta.ar1 ) 
      + f(tmean, model="rw2", hyper=theta.rw2 )
    fmly = c("binomial", "gamma")
    Z = Z.all
    
    R = inla( formula=fmla, family=fmly, 
             data=inla.stack.data(Z), 
             control.predictor=list(A=inla.stack.A(Z), compute=TRUE), 
             control.fixed=theta.beta.presence_abundance, 
             # control.family=theta.observations.presence_abundance,
             verbose=TRUE )
 

  
     
      
    # --------------------------
    # --------------------------
   

    # "GroupRho":: posterior marginal distribution of the temporal correlation
    (R$summary.hyperpar)
    round(R$summary.hy[3,], 5) 
    

    # random field parameters on user scale
    oo = inla.spde2.result(R, 'iAB', S0, do.transf=TRUE)
   
    coreltemp = pmatch( "GroupRho", names(R$marginals.hyper) )  # index of the temporal autocorrelation parameter

    plot(R$marginals.hyper[[coreltemp]], type='l', xlab=expression(rho), ylab='Density')  ## temporal autorcorrelation ("GroupRho") 
    # abline(v=rho, col=2)

    plot(oo$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[x]), ylab='Density')
    # abline(v=params[1], col=2)

    plot(oo$marginals.kappa[[1]], type='l', xlab=expression(kappa), ylab='Density')
    # abline(v=params[2], col=2)

    plot(oo$marginals.range.nominal[[1]], type='l', xlab='range nominal', ylab='Density')
    # abline(v=sqrt(8)/params[2], col=2)
   
    
    # indices for random field at data locations
    idat <- inla.stack.index( Z, 'iPA')$data

    # correlation between the the posterior mean and the response by
    cor( PA, R$summary.linear.predictor$mean[idat])

    [1] 0.761407





    # ----------------
    # prediction (of the latent field) for each time and visualize it .. first create projector from mesh to output

      pG = inla.mesh.projector( M0, xlim=p$corners$plon, ylim=p$corners$plat, dims=c(p$nplons, p$nplats)  )
      inside = inout( pG$lattice$loc, M0.domain$loc ) 
      
      xmean <- list()
      for (j in 1:nyrs) xmean[[j]] <- inla.mesh.project( pG, R$summary.random$i$mean[i$i.group==j])
      for (j in 1:nyrs) xmean[[j]][!inside] <- NA

      levelplot( xmean[[nyrs]], xlab='', ylab='', col.regions=topo.colors(100), scale=list(draw=FALSE), aspect="iso" )





    # Validation ...

      vdat <- data.frame(y=as.vector(y), w=ccov, time=rep(1:k, each=n), xcoo=rep(coords[,1], k), ycoo=rep(coords[,2], k))[-isel, ]
      Aval <- inla.spde.make.A(prmesh1, loc=cbind(vdat$xcoo, vdat$ycoo), group=vdat$time)
      stval <- inla.stack(tag='stval', data=list(y=NA), A=list(Aval,1), effects=list(iset, w=vdat$w))

      # Now, we just use a full data stack to fit the model
      stfull <- inla.stack(effdat, stval)

      vres <- inla(formulae, data=inla.stack.data(stfull),
          control.predictor=list(compute=TRUE, A=inla.stack.A(stfull)),
          control.family=list(initial=20, fixed=TRUE),
          control.inla=list(strategy='gaussian'),
          control.fixed=list(expand.factor.strategy='inla'))

      # We can look at fitted values for the validation data. We can plot the predicted versus
      # observed values to look at goodness of fit. First, we found the index for this data from
      # full stack data.


      ival <- inla.stack.index(stfull, 'stval')$data

      # We plot it with following commands and visualize at Figure 11.4.

      plot(vres$summary.fitted.values$mean[ival], vdat$y, asp=1, xlab='Posterior mean', ylab='Observed')
      abline(0:1, col=gray(.7))





---- Testing follows -- must tweak for correct variables, etc.

    
    pp = inla.tmarginal( function(x) 1/x, R$marginals.hyperpar [[1]] ) # sigma^2(e) .. nugget






    plot( R$marginals.fixed[[1]], type="l", xlab="Beta" )
    
       

    plot( R$marginals.fixed[["b0_PA"]], type="l", xlab="b0_PA" ) 
    plot( R$summary.random[["i"]][,c("ID", "mean")], type="l", xlab="" )
 
    
    str(R$marginals.hyperpar)
    plot( R$marginals.hyper[[1]], type="l", xlab="" )
    plot.default( inla.tmarginal( function(x) {1/exp(x)}, R$marginals.hyperpar[[1]]), xlab="", type="l")




