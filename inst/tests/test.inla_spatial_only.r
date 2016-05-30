
# try spatial only with inla

# Full spatial model with "Semicontinuous (Gamma, Bernoulli)" or delta-type response variable

# USEFUL FUNCTIONS (for INLA)

  II = function(x) {x}  # dummy function to return itself .. used below

  # see documentation here: http://www.r-inla.org/examples/tutorials/tutorial-on-option-scale-model
  inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors

# DATA
  # convert snow crab data into a spatial data frame
	p = snowcrab::initialise.local.environment( libs=c("mgcv", "lubridate", "lattice", "lattice", "grid", "fields", "parallel",
                         "sp", "INLA", "geostatsinla", "geostatsp", "raster"  ) )

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

  locs_all = as.matrix( set0[,c("plon", "plat")] )

  # boundary domain
  M0.domain = inla.nonconvex.hull( locs_all, convex=10, resolution=75 )


## DEBUG::
  testyear = 2011

  set0 = set0[ which(set0$yr == testyear ) , ]
  set0$plon = set0$plon
  set0$plat = set0$plat

  locs0  = as.matrix( set0[,c("plon", "plat")] )
  nData0 = nrow( set0 )


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
    i <- inla.spde.make.index('i', n.spde=S0$n.spde )
    iPA <- inla.spde.make.index('iPA', n.spde=S0$n.spde )
    iAB <- inla.spde.make.index('iAB', n.spde=S0$n.spde )

    # projection matrix A to translate from mesh nodes to data nodes
    A = inla.spde.make.A( mesh=M0, loc=locs0 )


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
      + f(iPA, model=S0, diagonal=1e-2)
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
      + f(iAB, model=S0, diagonal=1e-2 )
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
      + f( iAB, model=S0, diagonal=1e-4 )
      + f( iPA, copy='iAB', fixed=FALSE, diagonal=1e-4)
      + f( z, model="rw2", hyper=theta.rw2, diagonal=1e-4 )
      + f( tmean, model="rw2", hyper=theta.rw2, diagonal=1e-4 ) )
#    fmly = c("binomial", "gamma")
    fmly = c("binomial", "gaussian")
    Z = Z.all

    R = inla( formula=fmla, family=fmly,
             data=inla.stack.data(Z),
             control.predictor=list(A=inla.stack.A(Z), compute=FALSE ),
             control.fixed=theta.beta.presence_abundance,
             # control.family=theta.observations.presence_abundance,
             verbose=TRUE )

    summary(R)
    plot(R)



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


    (R$summary.hyperpar)


    # random field parameters on user scale
    oo = inla.spde2.result(R, 'iAB', S0, do.transf=TRUE)

    plot(oo$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[x]), ylab='Density')
    # abline(v=params[1], col=2)

    plot(oo$marginals.kappa[[1]], type='l', xlab=expression(kappa), ylab='Density')
    # abline(v=params[2], col=2)

    plot(oo$marginals.range.nominal[[1]], type='l', xlab='range nominal', ylab='Density')
    # abline(v=sqrt(8)/params[2], col=2)


    # indices for random field at data locations
    idat <- inla.stack.index( Z, 'iAB')$data

    # correlation between the the posterior mean and the response by
    cor( effdat$Y, R$summary.linear.predictor$mean[idat])

    [1] 0.761407





    # ----------------
    # prediction (of the latent field) for each time and visualize it .. first create projector from mesh to output

      pG = inla.mesh.projector( M0, xlim=p$corners$plon, ylim=p$corners$plat, dims=c(p$nplons, p$nplats)  )
      inside = inout( pG$lattice$loc, M0.domain$loc )

      xmean <- list()
      xmean <- inla.mesh.project( pG, R$summary.random$iAB$mean)
      xmean[!inside] <- NA

      levelplot( xmean, xlab='', ylab='', col.regions=topo.colors(100), scale=list(draw=FALSE), aspect="iso" )





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




