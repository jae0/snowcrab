
surplusproduction_model = function(  p, data, DS="stan" ) {


  if (DS=="stan" ) {

    if (0)  {
      library(rstan)
      current.year=2016
      p = bio.snowcrab::load.environment( year.assessment=current.year)
    }

    dir.create( p$surplusproduction_model$dir.output, recursive=T, showWarnings=F )
    fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"rdata", sep=".") )
    #fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"survey_final.rdata", sep=".") )

    sb = biomass.summary.db(p=p, DS="surplusproduction" )

    # ----------------------------------
    # these are already on log-scale
    sb$Kmu =  exp(sb$K.mu)
    sb$Ksd =  c(0.25, 0.25, 0.25) * sb$Kmu

    sb$rmu =  sb$r.mu
    sb$rsd =  c(0.25, 0.25, 0.25) * sb$rmu 

    sb$qmu =  sb$q.mu
    sb$qsd =  c(0.25, 0.25, 0.25) * sb$q.mu

    sb$missing = ifelse(is.finite(sb$IOA),0,1)
    sb$missing_n = colSums(sb$missing)
    sb$missing_ntot = sum(sb$missing_n)
    sb$IOA[  which(!is.finite(sb$IOA)) ] = 0  # reset NAs to 0 as stan does not take NAs


    surplus.stan = "
      data {
     
        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] qsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        matrix[N,U] CAT;
        matrix[N,U] IOA; 
        matrix[N,U] missing;
        int missing_n[U];
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        MN = M+N ;
        N1 = N+1;
      }

      parameters {
        vector <lower=eps>[U] K;
        vector <lower=eps,upper=3>[U] r;
        vector <lower=eps,upper=2>[U] q;
        vector <lower=eps,upper=2>[U] qs;
        vector <lower=eps,upper=(1-eps)>[U] bosd;  // observation error
        vector <lower=eps,upper=(1-eps)>[U] bpsd;  // process error
        vector <lower=eps,upper=(1-eps)>[U] b0;
        vector <lower=eps>[missing_ntot] IOAmissing;
        matrix <lower=eps>[M+N,U] bm;
      }
      
      transformed parameters {
        matrix[N,U] Y;  // index of abundance 
        matrix[N,U] Ymu;  // collator used to force positive values for lognormal
        matrix[N,U] logYmu;

        matrix[MN,U] bmmu; // collator used to force positive values for lognormal
        matrix[MN,U] logbmmu; 
        matrix[MN,U] rem;  // observed catch

        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
              Y[i,j] = IOA[i,j];
              if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }

        // -------------------  
        // removals (catch) observation model, standardized to K (assuming no errors in observation of catch!)
        for (j in 1:U) {
          rem[1:N,j] =  CAT[1:N,j]/K[j] ;
          rem[(N+1):MN,j] =  er*bm[ N:(MN-1),j] ;  // forecasts
        }

        // -------------------  
        // observation model calcs and contraints: 
        // Ymu = 'surveyed/observed' residual biomass at time of survey (Bsurveyed)
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t) 
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1) 
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line: 
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        // assume similar to a spring fishery but no need for separate q's
        //    Btot(t) = Bsurveyed(t)+ removals(t-1)
        //    NOTE: year designation in 4X is for the terminal year: ie. 2001-2002 => 2002

        for (j in 1:2) {  
          Ymu[1,j]        = qs[j] * bm[1,j] - rem[1,j] ; // starting year approximation
          Ymu[2:(ty-1),j] = qs[j] * bm[2:(ty-1),j] - rem[1:(ty-2),j] ; //spring surveys
          Ymu[ty,j]       = q[j]  * bm[ty,j] - (rem[(ty-1),j] + rem[ty,j] )/2.0  ; // transition year .. approximation
          Ymu[(ty+1):N,j] = q[j]  * bm[(ty+1):N,j] - rem[(ty+1):N,j] ;   // fall surveys    
        }
        {
          int k;
          k=3;
          Ymu[1,k]        = qs[k] * bm[1,k]   - rem[1,k] ; // starting year approximation
          Ymu[2:(ty-1),k] = qs[k] * bm[2:(ty-1),k] - rem[1:(ty-2),k]; 
          Ymu[ty:N,k]     = q[k]  * bm[ty:N,k] - rem[(ty-1):(N-1),k]; 
        }
        
        for (j in 1:U) {
          for (i in 1:N) {
            Ymu[i,j] = K[j] * fmax( Ymu[i,j], eps); // force positive value
          }
        }


        // -------------------  
        // process model calcs and constraints
        for (j in 1:U) {
          bmmu[1,j] = b0[j] ; // biomass at first year
          for (i in 2:MN) {
            bmmu[i,j] = bm[i-1,j] * ( 1.0 + r[j]*(1-bm[i-1,j]) ) - rem[i-1,j] ;
          }
        }
        for (j in 1:U) {
          for (i in 1:MN) {
            bmmu[i,j] = fmax(bmmu[i,j], eps);  // force positive value
          }
        }

        logYmu = log( Ymu ) ;
        logbmmu = log( bmmu );

      }

      model {

        // -------------------  
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        q ~ normal( qmu, qsd )  ;
        qs ~ normal( qmu, qsd )  ;
        b0 ~ beta( 8, 2 ) ; // starting b prior to first catch event 
        bosd ~ beta( 2, 8 ) ;
        bpsd ~ beta( 2, 8 ) ;



        // -------------------  
        // biomass observation model 
        for (j in 1:U) {  
          Y[1:N,j] ~ lognormal( logYmu[1:N,j], bosd[j] ) ;  
        }

        // -------------------  
        // biomass process model 
        for (j in 1:U) {
          bm[1:MN,j] ~ lognormal( logbmmu[1:MN,j], bpsd[j] ) ;     
        }
      
      }

      generated quantities {
        // matrix[MN,U] pd;
        // vector[U] MSY;
        // vector[U] BMSY;
        // vector[U] FMSY;
        // matrix[MN,U] B;
        // matrix[MN,U] P;
        // matrix[MN,U] C;
        
        // matrix[MN,U] F;
        // matrix[M,U] TAC;


        // -------------------  
        // annual production
        // for(j in 1:U) {
        //   pd[1,j] = bm[2,j]- bm[1,j] + rem[1,j] ; // approximation
        //   for (i in 2:N ){
        //     pd[i,j] = (bm[i+1,j]- bm[i-1,j])/2 + rem[i,j] ; // linear interpolation cancels out the bm[i,j] term
        //   }
        //   for(i in N1:(MN-1)) {
        //     pd[i,j] = (bm[i+1,j]- bm[i-1,j])/2 + er * bm[i-1,j] ;  // linear interpolation cancels out the bm[i,j] term
        //   }
        //   pd[MN,j] = (bm[MN,j]- bm[(MN-1),j]) + er * bm[(MN-1),j]  ; // approximation
        // }

        // -------------------  
        // fishing mortality
        // force first year estimate assuming catches in year 0 to be similar to year 1 

        // for (j in 1:U) {
        //   F[1,j] =  1.0 - rem[1,j] / bm[1,j] ;
        //   for (i in 2:MN) {
        //     F[i,j] =  1.0 - er * bm[i-1,j] / bm[i,j]  ;
        //   }
        // }
        // for (j in 1:U) {
        //   for (i in 1:MN) {
        //     F[i,j] =  -log( fmax( F[i,j], eps) )  ;
        //   }
        // }
     
        // -------------------  
        // parameter estimates for output
        // for(j in 1:U) {
        //   MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
        //   BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
        //   FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
      //    BX2MSY[j] = 1.0 - step( bm[N1,j]-0.25 ) ; // test if bm >= 1/2 bmY
      //    Bdrop[j]  = 1.0 - step( bm[N1,j]-bm[N,j] ) ; // test if bm(t) >= bm(t-1) 
      //    Fcrash[j] = 4.0 * MSY[j] / exp(K[j]) ; // fishing mortality at which the stock will crash
        // }

        // recaled estimates
        // for(j in 1:U) {
        //   for(i in 1:MN) {
        //     B[i,j] = (bm[i,j] - rem[i,j]) * K[j] ;
        //     P[i,j] = pd[i,j]*K[j] ;
        //     C[i,j] = rem[i,j]*K[j] ;
        //   }
        //   for(i in 1:M) {
        //     TAC[i,j] = rem[N+i,j]*K[j] ;
        //   }
        // }

      }

    "


    stanmodel = rstan::stan_model( model_code=surplus.stan )

    f = sampling(stanmodel, data=sb, chains=5, iter=5000, warmup = 1000,
      control = list(adapt_delta = 0.85, max_treedepth=12) )
          # warmup = 200,          # number of warmup iterations per chain
          # control = list(adapt_delta = 0.9),
          # # refresh = 500,          # show progress every 'refresh' iterations
          # iter = 1000,            # total number of iterations per chain
          # chains = 5,             # number of Markov chains
          # cores = 5              # number of cores (using 2 just for the vignette)


    plot(f)
    print(f)
    traceplot(f)

    # extract samples
    e = rstan::extract(f, permuted = TRUE) # return a list of arrays
    m2 = as.array(f)

    traceplot(f, pars=c("K"))
    pred=rstan::extract(f)

    est=colMeans(pred)
  
    prob=apply(pred,2,function(x) I(length(x[x>0.10])/length(x) > 0.8)*1)

  }


  if (DS=="jags") {

    require(rjags)
    rjags::load.module("dic")
    rjags::load.module("glm")



    dir.create( p$surplusproduction_model$dir.output, recursive=T, showWarnings=F )
    fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"rdata", sep=".") )
    #fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"survey_final.rdata", sep=".") )

    sb = biomass.summary.db(p=p, DS="surplusproduction" )

    sb$tomonitor = c( "r", "K", "q", "qs", "r.mu", "r.sd", "b","bp.sd", "bo.sd","b0", "b0.sd", "rem", "rem.sd", "rem.mu","REM", "MSY", "BMSY", "FMSY", "Fcrash", "Bdrop", "BX2MSY", "F", "TAC",  "C", "P", "B" )

    sb$jagsmodelname = "biomassdynamic_nonhyper_2016.bugs"


    m = jags.model( file=fishery.model.jags ( DS=sb$jagsmodelname ), data=sb, n.chains=n.chains, n.adapt=n.adapt ) # recruitment + spring/summer q's + all observed CVs
    tomonitor = intersect( variable.names (m), sb$tomonitor )
    y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.final, thin=n.thin) # sample from posterior
    dic.samples(m, n.iter=n.iter ) # pDIC

    graphics.off() ; x11(); layout( matrix(c(1,2,3), 3, 1 )); par(mar = c(5, 4, 0, 2))
    for( i in 1:3) hist(y$cv.r[i,,], "fd")

    # convergence testing -- by 1000 to 1500 convergence observed by Gelman shrink factor diagnostic
    y = jags.samples(m, variable.names=tomonitor, n.iter=6000, thin=1 )

    gelman.plot(y[["r"]])
    gelman.plot(y[["K"]])
    gelman.plot(y[["q"]])  # about 6-8000 runs required to converge
    gelman.plot(y[["r.sd"]])
    gelman.plot(y[["K.sd"]])
    gelman.plot(y[["bo.sd"]])
    gelman.plot(y[["bp.p"]])
    geweke.plot(y[["r"]])

    # update if not yet converged
    #  update(m, n.iter=n.iter ) # above seems enough for convergence but a few more to be sure

    # determine autocorrelation thinning
    y = coda.samples(m, variable.names=c("K", "r", "q"), n.iter=20000, thin=10) # sample from posterior
    autocorr.plot(y)   # about 10 to 20 required
    # plot(y, ask=T)
    # autocorr(y, lags = c(0, 1, 5, 10, 50), relative=TRUE)

    # final sampling from the posteriors
    #  y = jags.samples(m, variable.names=tomonitor, n.iter=10000, thin=20) # sample from posterior
    y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.final, thin=n.thin) # sample from posterior
    save(y, file=fnres, compress=T)
    # load( fnres )



    ## 
    b1=apply( y$B[,1,,], 1,quantile , probs=c(0.025,0.5,0.975), na.rm=T  )
    f1=apply( y$F[,1,,], 1,quantile , probs=c(0.5), na.rm=T  )
    NENS=data.frame(t(b1),F=f1,row.names=1999:2019)
    NENS$U=1-exp(-NENS$F)
    b2=apply( y$B[,2,,], 1,quantile , probs=c(0.025,0.5,0.975), na.rm=T  )
    f2=apply( y$F[,2,,], 1,quantile , probs=c(0.5), na.rm=T  )
    SENS=data.frame(t(b2),F=f2,row.names=1999:2019)
    SENS$U=1-exp(-SENS$F)

    b3=apply( y$B[,3,,], 1,quantile , probs=c(0.025,0.5,0.975), na.rm=T  )
    f3=apply( y$F[,3,,], 1,quantile , probs=c(0.5), na.rm=T  )
    CFA4X=data.frame(t(b3),F=f3,row.names=1999:2019)
    CFA4X$U=1-exp(-CFA4X$F)

    NENS
    SENS
    CFA4X
    # Figures
    graphics.off()


    # frequency density of key parameters
    figure.bugs( "K", y=y, sb=sb, fn=file.path(dir.output, "K.density.png" ) )
    figure.bugs( "r", y=y, sb=sb, fn=file.path(dir.output, "r.density.png" ) )
    figure.bugs( "q", y=y, sb=sb, fn=file.path(dir.output, "q.density.png" ) ,xrange=c(0,2))
    figure.bugs( "FMSY", y=y, sb=sb, fn=file.path(dir.output, "FMSY.density.png" ) )
    figure.bugs( "bo.sd", y=y, sb=sb, fn=file.path(dir.output, "bo.sd.density.png" ) )
    figure.bugs( "bp.sd", y=y, sb=sb, fn=file.path(dir.output, "bp.sd.density.png" ) )

    # timeseries
    figure.bugs( type="timeseries", vname="biomass", y=y, sb=sb, fn=file.path(dir.output, "biomass.timeseries.png" ), save.plot=T )
    figure.bugs( type="timeseries", vname="fishingmortality", y=y, sb=sb, fn=file.path(dir.output, "fishingmortality.timeseries.png" ) )

    # Harvest control rules
    figure.bugs( type="hcr", vname="default", y=y, sb=sb, fn=file.path(dir.output, "hcr.default.png" ), save.plot=T  )
    figure.bugs( type="hcr", vname="simple", y=y, sb=sb, fn=file.path(dir.output, "hcr.simple.png" ) )

    # diagnostics
    figure.bugs( type="diagnostic.production", y=y, sb=sb, fn=file.path(dir.output, "diagnostic.production.png" ) )
    figure.bugs( type="diagnostic.errors", y=y, sb=sb, fn=file.path(dir.output, "diagnostic.errors.png" ) )
    figure.bugs( type="diagnostic.phase", y=y, sb=sb, fn=file.path(dir.output, "diagnostic.phase.png" ) )


    ndata = sb$N

    # densities of biomass estimates for the year.assessment
      for (i in 1:3) plot(density(y$B[ndata,i,,] ), main="")
      qs = apply( y$B[ndata,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
      qs

      # densities of biomass estimates for the previous year
      for (i in 1:3) plot(density(y$B[ndata-1,i,,] ), main="")
      qs = apply( y$B[ndata-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
      qs

      # densities of F in assessment year
      for (i in 1:3) plot(density( y$F[ndata,i,,] ), xlim=c(0.05, 0.5), main="")
      qs = apply( y$F[ndata,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
      qs
      qs = apply( y$F[ndata,,,], 1, mean )

      # densities of F in previous year
      for (i in 1:3) plot(density( y$F[ndata-1,i,,] ), main="")
      qs = apply( y$F[ndata-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
      qs

      # F for table ---
      summary(y$F, median)


    debug = F
    if (debug) {

      graphics.off()

      x11()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(5, 4, 0, 2))

      for( i in 1:3) hist(y$r[i,,], "fd")
      for( i in 1:3) hist(y$q.sd[i,,], "fd")
      for( i in 1:3) hist(y$K.sd[i,,], "fd")
      for( i in 1:3) hist(y$b0.sd[i,,], "fd")
      for( i in 1:3) hist(y$r.sd[i,,], "fd")
      for( i in 1:3) hist(y$b0[i,,], "fd")

    }



  }


  if (DS=="LaplacesDemon") {
    require(LaplacesDemonCpp)
      
    if (!exists("current.year")) current.year=lubridate::year(Sys.Date())
    p = bio.snowcrab::load.environment( year.assessment=current.year)


    dir.create( p$surplusproduction_model$dir.output, recursive=T, showWarnings=F )
    fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"rdata", sep=".") )
    #fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"survey_final.rdata", sep=".") )

    sb = biomass.summary.db(p=p, DS="surplusproduction" )

    debug.region="cfa4x"
    debug.region="cfasouth"
    debug.region="cfanorth" 

    # single region test
    K0est =list()
    K0est[["cfanorth"]] = 5
    K0est[["cfasouth"]] = 65
    K0est[["cfa4x"]] = 5
        
    yrs = as.numeric(rownames( res$B))
    nforecasts = 5
    

    Data = list(
      Ndata = length(yrs),
      Nforecasts = nforecasts, # no years for projections
      O = res$B[[debug.region]], # observed index of abundance
      log_O0 = mean( log(res$B[[debug.region]]), na.rm=TRUE ),
      Omax = max( res$B[[debug.region]] *1.5 , na.rm=TRUE),
      removals = res$L[[debug.region]] , # removalsches  , assume 20% handling mortality and illegal landings
      log_R0 = log(mean(res$L[[debug.region]], na.rm=TRUE )/K0est[[debug.region]]),
      er = 0.2,  # target exploitation rate
      ty = which(yrs==2004) ,  # index of the transition year (2004) between spring and fall surveys
      log_r0= log(1),
      log_K0= log(K0est[[debug.region]]),
      log_q0= log(mean(res$B[[debug.region]], na.rm=TRUE)/K0est[[debug.region]]),
      S0= 0.6, # normalised
      cv = 0.4,
      smax =1.25,
      eps = 1e-6, 
      eps_1 = 1 - 1e-6
    )



    # set up the model
  # set up model for a simple surplus production model 
  # to be solved by the Rlibrary: LaplacesDemon (or alternately via penalized Maximum Likelihood) 

  

      # ----------------------------------
      # Identify location and number of missing values -- prediction locations are treated the same way

      # compute a few things here that are constant in the model
      Data$N = Data$Ndata + Data$Nforecasts # no years with data + projections
      Data$q0 = exp( Data$log_q0) 
      Data$r0 = exp( Data$log_r0) 
      Data$K0 = exp( Data$log_K0) 
      Data$R0 = exp(Data$log_R0)
      Data$O0 = exp(Data$log_O0)
      Data$O_range = range( Data$O, na.rm=TRUE)

      # ----------------------------------

      Data$Missing = list(
        O = list( 
          n=length(intersect( 1:Data$Ndata, which( !is.finite(Data$O) ) )), 
          idx=intersect( 1:Data$Ndata, which( !is.finite(Data$O) ) ) 
        ), # observations of survey index
        removals = list( 
          n=length(intersect( 1:Data$Ndata, which( !is.finite(Data$removals) ) ) ), 
          idx=intersect( 1:Data$Ndata, which( !is.finite(Data$removals) ) )  
        )  # observations of removals
      )

      # ----------------------------------

      if ( Data$Nforecasts > 0 ) {
        # extend data for forecasts
        Data$O = c( Data$O, rep(NA, Data$Nforecasts) )
        Data$removals = c( Data$removals, rep(NA, Data$Nforecasts) )
        Data$Forecasts = list( 
          n = Data$Nforecasts,
          idx = Data$Ndata + 1:Data$Nforecasts,
          idx_1 = (Data$Ndata + 1:Data$Nforecasts) -1 
        )
      }

      # ----------------------------------
      # misc indexing
      Data$idx = list(
        t1 = 1,
        t2 = 2:(Data$ty-1),
        t2_1 = (2:(Data$ty-1) ) - 1,
        t3 =  Data$ty,
        t3_1 =  Data$ty-1,
        t4 =  (Data$ty+1):Data$N,
        icurrent  = 2:Data$N,
        iprevious = 1:(Data$N-1),
        idata = 1:Data$Ndata
      )

      # ----------------------------------
      # paramater names and initial values
      Data$parm.names = list(
        r = Data$r0,
        K = Data$K0,
        q = Data$q0,
        r_sd = Data$r0,
        K_sd = Data$K0,
        q_sd = Data$q0,
        S0_sd = Data$S0,
        S_sd = Data$S0,
        O_sd = Data$O0,
        R_sd = Data$R0,
        S = rep( Data$S0, Data$N),
        S0 = Data$S0
      ) 
      if (Data$Missing$removals$n > 0) Data$parm.names$removals_imp = rep(0, Data$Missing$removals$n)
      if (Data$Missing$O$n > 0) Data$parm.names$O_imp = rep(0, Data$Missing$O$n)
      Data$parm.names = as.parm.names( Data$parm.names )

      # ----------------------------------
      # index position of paramaters
      Data$pos = list(
        r=grep("\\<r\\>", Data$parm.names),
        K=grep("\\<K\\>", Data$parm.names),
        q=grep("\\<q\\>", Data$parm.names),
        r_sd = grep("\\<r_sd\\>", Data$parm.names), 
        K_sd = grep("\\<K_sd\\>", Data$parm.names),
        q_sd = grep("\\<q_sd\\>", Data$parm.names),
        S0_sd = grep("\\<S0_sd\\>", Data$parm.names),
        S_sd=grep("\\<S_sd\\>", Data$parm.names),
        O_sd=grep("\\<O_sd\\>", Data$parm.names),
        R_sd =grep("\\<R_sd\\>", Data$parm.names),
        S=grep("\\<S\\>", Data$parm.names),
        S0=grep("\\<S0\\>", Data$parm.names)
      )
      if (Data$Missing$removals$n > 0) Data$pos$removals_imp = grep("\\<removals_imp\\>", Data$parm.names)
      if (Data$Missing$O$n > 0) Data$pos$O_imp = grep("\\<O_imp\\>", Data$parm.names)

      # ----------------------------------
      # monitoring nodes
      Data$mon.names = c("LP", "r", "K", "q" )
      # Data$mon.names = c("LP", "r", "K", "q", paste0("S",1:Data$N), paste0("AR",1:(Data$N-1) ) )


      # ----------------------------------
      # Parameter Generating Function
      # these parameters are operate on the log-scale to force positive values ...
      Data$param$means = c( 
        r=Data$r0, K=Data$K0, q=Data$q0, 
        r_sd=Data$r0/2, K_sd=Data$K0/2,  q_sd=Data$q0/2, 
        S0_sd=Data$S0/2, S_sd=0.2, O_sd=Data$O0/2, R_sd=Data$R0/2,
        S=rep( 0.5, Data$N), S0=Data$S0
      )
      if (Data$Missing$removals$n > 0) {
        Data$param$means = c( Data$param$means, removals_imp=rep(Data$R0, Data$Missing$removals$n ) )
      }
      if (Data$Missing$O$n > 0) {
        Data$param$means = c( Data$param$means, O_imp=rep(Data$O0, Data$Missing$O$n) )
      }
      Data$param$sd = Data$param$means / 2 
      Data$param$n = length(Data$param$means)

      Data$PGF = function(Data) {
        out = LaplacesDemonCpp::rnorm(Data$param$n, Data$param$means, Data$param$sd )
        i = which( out < Data$eps )
        if (length(i)>0) out[i] = -out[i] # reflect onto positive side
        return( out  )
      }
      Data$PGF = compiler::cmpfun(Data$PGF)


      # ----------------------------------
      # define the model that generates the loglikelihoods (LL) and the log-posteriors (LP)
      Data$Model = function(parm, Data) {
        
        Spred = Opred = rep(0, Data$N )   # initialize a few storage vectors
          # NOTE: for lognormal: cv = sqrt(exp(sigma^2) - 1); 
        # or sigma = sqrt(log(cv^2+ 1) ) ==> sigma = sqrt( log(0.25^2 + 1)) = 0.246 ~ cv -- i.e. cv ~ sd
        i = which( parm < Data$eps )
        if (length(i)>0) parm[i] = -parm[i] # reflect onto positive side

        q =  parm[Data$pos$q]
        r =  parm[Data$pos$r]
        K =  parm[Data$pos$K]
        S =  parm[Data$pos$S]
        Spred[1] =  parm[Data$pos$S0]
       
       # SD params stay on log-scale
        q_sd = parm[Data$pos$q_sd]
        r_sd = parm[Data$pos$r_sd]
        K_sd = parm[Data$pos$K_sd]
        S0_sd = parm[Data$pos$S0_sd]
        O_sd = parm[Data$pos$O_sd]
        S_sd = parm[Data$pos$S_sd]
        R_sd = parm[Data$pos$R_sd]
        
        llkimpute = 0
        llkprior = c()
        llkprior[Data$parm.names] = 0
        llkprior[Data$pos$q] = dnorm( log(parm[Data$pos$q]), Data$log_q0, q_sd, log=TRUE ) ;
        llkprior[Data$pos$r] = dnorm( log(parm[Data$pos$r]), Data$log_r0, r_sd, log=TRUE ) ;
        llkprior[Data$pos$K] = dnorm( log(parm[Data$pos$K]), Data$log_K0, K_sd, log=TRUE ) ;
        llkprior[Data$pos$S0] = dnorm( log(parm[Data$pos$S0]), log(Data$S0), S0_sd, log=TRUE ) ;
        llkprior[Data$pos$q_sd] = dgamma( q_sd, shape=Data$q0, scale=100, log=TRUE );
        llkprior[Data$pos$r_sd] = dgamma( r_sd, shape=Data$r0, scale=100, log=TRUE );
        llkprior[Data$pos$K_sd] = dgamma( K_sd, shape=Data$K0, scale=100, log=TRUE );
        llkprior[Data$pos$S0_sd] = dgamma( S0_sd, shape=Data$S0, scale=100, log=TRUE );
        llkprior[Data$pos$O_sd] = dgamma( O_sd, shape=Data$O0, scale=100, log=TRUE );
        llkprior[Data$pos$S_sd] = dgamma( S_sd, shape=0.5, scale=100, log=TRUE );
        llkprior[Data$pos$R_sd] = dgamma( R_sd, shape=Data$R0, scale=100, log=TRUE );

        R = Data$removals/K ;# make sure it is producing sensible values:
        # impute missing data 
        if ( Data$Missing$removals$n > 0 ) {
          parm[Data$pos$removals_imp] = R[Data$Missing$removals$idx ] = LaplacesDemonCpp::interval( parm[Data$pos$removals_imp], Data$eps, Data$smax)  
          llkimpute = llkimpute + sum( dnorm( log(R[Data$Missing$removals$idx ]), Data$log_R0, R_sd, log=TRUE ) );
        }

        if ( Data$Missing$O$n > 0 ) {
          parm[Data$pos$O_imp] = Data$O[Data$Missing$O$idx ] = LaplacesDemonCpp::interval( parm[Data$pos$O_imp], Data$eps, Data$K0)  
          llkimpute = llkimpute + sum(dnorm( log(Data$O[Data$Missing$O$idx ]), Data$log_O0, O_sd, log=TRUE ) ) ;    
        }

        if ( Data$Forecasts$n > 0 ) {
          R[Data$Forecasts$idx] = S[Data$Forecasts$idx_1] * Data$er 
          llkimpute = llkimpute + sum(dnorm( log(R[Data$Forecasts$idx]), Data$log_R0, R_sd, log=TRUE )) ;
          Data$O[Data$Forecasts$idx] = runif( Data$Forecasts$n, Data$O_range[1], Data$O_range[2]  )
          # likelihood for Data$O is below
        }

        # process model -- simple logistic
        Spred[ Data$idx$icurrent] = S[Data$idx$iprevious] * 1.0 + r*{1-S[Data$idx$iprevious]} - R[Data$idx$iprevious] ;  # 
        Spred = LaplacesDemonCpp::interval_random( Spred, Data$eps, Data$smax, 0.01 ) 
        llkprior[Data$pos$S] = dnorm( log(parm[Data$pos$S]), log(Spred), S_sd, log=TRUE ) 

        # observation model
        Opred[Data$idx$t1] = S[Data$idx$t1] - R[Data$idx$t1] ; # first year approximation
        Opred[Data$idx$t2] = S[Data$idx$t2] - R[Data$idx$t2_1] ;
        Opred[Data$idx$t3] = S[Data$idx$t3] - (R[Data$idx$t3_1] + R[Data$idx$t3])/2 ; # transition year from Spring to Autumn survey
        Opred[Data$idx$t4] = S[Data$idx$t4] - R[Data$idx$t4] ;
        Opred = K*q*Opred
        Opred = LaplacesDemonCpp::interval_random( Opred, Data$eps, Inf, 0.01 )
        llkdata = dnorm( log(Data$O), log(Opred), O_sd, log=TRUE )
     
        # additional computed variables of interest 
        ER = R / S ;
        ER = LaplacesDemonCpp::interval_random( ER, Data$eps, Data$eps_1, 0.01 ) 
        B = S*K
        C = R*K
        F = -log( 1 - ER) ; # fishing mortality
        AR = S[ Data$idx$icurrent] / S[Data$idx$iprevious] ;  # "autocorelation" component

        LL = sum( llkdata )  # log likelihood (of the data)
        LP = LL + sum( llkprior ) + sum(llkimpute) # log posterior
        out = list( LP=LP, Dev=-2*LL, Monitor=c(LP, r, K, q), yhat=S*K, parm=parm )
        return( out  )
      }

      Data$Model.ML  = compiler::cmpfun( function(...) (Data$Model(...)$Dev / 2) )  # i.e. - log likelihood
      Data$Model.PML = compiler::cmpfun( function(...) (- Data$Model(...)$LP) ) #i.e., - log posterior 
      Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed

      print (Data$Model( parm=Data$PGF(Data), Data=Data ) ) # test to see if return values are sensible
     


    # 2. maximum likelihood solution
    f.ml = optim( par=sb$PGF(sb), fn=sb$Model.ML, Data=sb, method="BFGS", control=list(maxit=5000, trace=0), hessian=TRUE  )
    names(f.ml$par ) = sb$parm.names
    #print(sqrt( diag( solve(f.ml$hessian) )) ) # assymptotic standard errors
    (f.ml$par)

    # 3. penalized maximum likelihood .. better but still a little unstable depending on algorithm
    f.pml = optim( par=sb$PGF(sb), fn=sb$Model.PML, Data=sb, method="BFGS", control=list(maxit=5000, trace=0), hessian=TRUE )
    names(f.pml$par ) = sb$parm.names
    (f.pml$par)

    #print(sqrt( diag( solve(f.pml$hessian) )) ) # assymptotic standard errors



    f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=1000, Status=100, Thinning=10) # 

    f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=10) # 

    f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=as.initial.values(f), Iterations=5000, Status=10, Thinning=20, Method="NUTS") #


    Initial.Values <- as.initial.values(f)
    f <- LaplacesDemon(sb$Model, Data=sb, Initial.Values,
         Covar=f$Covar, Iterations=20000, Status=1000, Thinning=1000,
         Algorithm="CHARM", Specs=NULL)


    f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=as.initial.values(f), Iterations=5000, Status=100, Thinning=25, Covar=f$Covar)  



    f = LaplaceApproximation(sb$Model, Data=sb, parm=sb$PGF(sb), Method="Roptim", method="BFGS", Stop.Tolerance=1e-9, Iterations = 10000  )
    f

    Consort(f)
    plot(f, Data=sb)

    PosteriorChecks(f)


    f0 = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=1000, Status=100, Thinning=1)


    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="HAR", Iterations=1000  ) 

    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="SR1", Iterations=1000  ) 

    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="TR", Iterations=10000  ) 


    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="BFGS", Iterations=1000  ) 


    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="PSO", Iterations=10000  ) 


    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f0), Method="SPG", Iterations=10000  ) 




    # 4. quick solution: acts as "burn-in" .. do a few times in case solution is unstable
    f = LaplaceApproximation(sb$Model, Data=sb, parm=sb$PGF(sb), Iterations=100 ) 
    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Iterations=1000 ) 
    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="BFGS", Stop.Tolerance=1e-8, Iterations=1000  ) 


    f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=5000, Status=10, Thinning=10, Algorithm="NUTS", Covar=f$Covar, Specs=list(A=100, delta=0.6, epsilon=NULL, Lmax=Inf))  # A=burnin, delta=target acceptance rate

    f = LaplacesDemon.hpc(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=10, Status=1, Thinning=1, Algorithm="NUTS", Covar=f$Covar, Specs=list(A=100, delta=0.6, epsilon=NULL, Lmax=Inf))  # A=burnin, delta=target acceptance rate

    f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=1000, Status=100, Thinning=1)

    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="Roptim", Stop.Tolerance=1e-8, Iterations = 5000, method="Nelder-Mead", control=list(maxit=5, reltol=1e-9 ) )
    f

    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="Roptim", Stop.Tolerance=1e-8, Iterations = 1000, method="L-BFGS-B", control=list(maxit=10, reltol=1e-10 ) )



    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f)  ) # fast spin up of paramters
    f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="BFGS", Stop.Tolerance=1e-8  ) 


    # NOTES: NM is slow
    # algorithms = c( "TR", "LM", "CG", "NR", "NM", "SPG", "LBFGS", "BFGS", "PSO", "SR1", "HAR", "DFP", "BHHH", "HJ", "Rprop" ) # available
    # algorithms = c( "CG", "NM", "SPG", "LBFGS", "BFGS", "PSO", "SR1", "HAR", "DFP", "Rprop" ) # reliably working
    f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=500, Status=100, Thinning=1)

    algorithms = c( "CG", "LBFGS", "BFGS", "PSO", "SR1", "HAR", "DFP", "Rprop" ) # reliably working
    for (a in algorithms) {
      print(a)
      ft = try( LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method=a, Stop.Tolerance=1e-9, Iterations=5000 ) )
      if (! class(ft) %in% "try-error" ) if (ft$Converged) if( ft$LP.Final > f$LP.Final )  {
        f = ft
        f$Method = a 
      }
    }
    str(f)
    plot(f, Data=sb)

    # MCMC ... look at ACF
    f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=as.initial.values(f), Covar=f$Covar, 
      Iterations=1000, Status=100, Thinning=1, Algorithm="CHARM" )
    Consort(f)
    plot(f, Data=sb)

     PosteriorChecks(f)
         #caterpillar.plot(f, Parms="beta")
         #plot(f, Data, PDF=FALSE)
         #Pred <- predict(f, Model, Data, CPUs=1)

      #summary(Pred, Discrep="Chi-Square")
         #plot(Pred, Style="Covariates", Data=Data)
         #plot(Pred, Style="Density", Rows=1:9)
         #plot(Pred, Style="Fitted")
         #plot(Pred, Style="Jarque-Bera")
         #plot(Pred, Style="Predictive Quantiles")
         #plot(Pred, Style="Residual Density")
         #plot(Pred, Style="Residuals")
         #Levene.Test(Pred)
         #Importance(f, Model, Data, Discrep="Chi-Square")

    # MCMC: run with appropriate thinning and other options:
    f = LaplacesDemon(sb$Model, Data=sb, as.initial.values(f),
      Covar=NULL, Iterations=10000, Status=1000, Thinning=1000, Algorithm="CHARM", Specs=NULL)

    f = LaplacesDemon(sb$Model, Data=sb, as.initial.values(f),
      Covar=f$Covar , Iterations=50000, Status=10204, Thinning=1000, Algorithm="CHARM", Specs=list(alpha.star=0.44))


    f = LaplacesDemon(sb$Model, Data=sb, 
      Initial.Values=as.initial.values(f), Covar=f$Covar, Iterations=10000, Status=1000, Thinning=100)

    f = LaplacesDemon(sb$Model, Data=sb, 
      Initial.Values=as.initial.values(f), Covar=f$Covar, Iterations=30000, Status=1000, Thinning=35, 
      Algorithm="AMWG", Specs=list(B=NULL, n=1000, Periodicity=35 ) )

    f <- LaplacesDemon(sb$Model, Data=sb,
      Initial.Values=as.initial.values(f), Covar=f$Covar, Iterations=10000, Status=1000, Thinning=105, 
      Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))

    # medium speed .. increase Iterations til convergence
    f = VariationalBayes(sb$Model, Data=sb,  parm=as.initial.values(f), Iterations=100,  Samples=10, CPUs=5 )
    f = VariationalBayes(sb$Model, Data=sb,  parm=as.initial.values(f), Iterations=5000,  Samples=100, CPUs=5 ,Stop.Tolerance=1e-9)

    # slow
    f = IterativeQuadrature(sb$Model, Data=sb, parm=as.initial.values(f), Iterations=10, Algorithm="AGH",
     Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )



  }


  if (DS=="tmb") {


    require(bio.base)
      
    if (!exists("current.year")) current.year=lubridate::year(Sys.Date())
    p = bio.snowcrab::load.environment( year.assessment=current.year)


    dir.create( p$surplusproduction_model$dir.output, recursive=T, showWarnings=F )
    fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"rdata", sep=".") )
    #fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"survey_final.rdata", sep=".") )

    sb = biomass.summary.db(p=p, DS="surplusproduction" )


    debug.region="cfa4x"
    debug.region="cfanorth" 
    debug.region="cfasouth"


    # Data
    data=list(N=sb$Ndata, IOA=sb$O, CAT=sb$removals);
    # Parameters: initial values
    parameters=list(log_sigmap=log(0.2),
                    log_sigmao=log(0.2),
                    log_Q=log(0.1),
                    log_r=log(0.8),
                    log_K=log(70),
                    log_P=rep(log(0.6),sb$Ndata))

    #parameters=list(log_tau=rep(log(0.1),sb$U),
    #                log_sigma=rep(log(0.1),sb$U),
    #                log_Q=rep(log(0.1),sb$U),
    #                log_r=rep(log(0.8),sb$U),
    #                log_K=log(c(6,70,2)),
    #                log_P=matrix(rep(log(0.6),sb$N*sb$U),sb$N,sb$U))
    #

    ############################################################################
    # Compile model cpp file in TMB
    ############################################################################

    tmb.dir = file.path(project.codedirectory('bio.snowcrab'),"inst","tmb")
    setwd(tmb.dir)

    library(TMB);
    compile("biomassdynamic.cpp")
    dyn.load(dynlib("biomassdynamic"))

    ###########################################################################
    # Estimation
    ###########################################################################

    # The objective function
    obj <- MakeADFun(data, parameters,  random="log_P", DLL="biomassdynamic")  

    # Optimize the objective function
    opt <- nlminb(obj$par,obj$fn,obj$gr)

    # Report estimates and standard errors
    rep <- sdreport(obj)

    MLE <- list(year=1:data$N,B=rep$value[1:data$N], B.sd=rep$sd[1:data$N],predI=rep$value[1:data$N+data$N], predI.sd=rep$sd[1:data$N+data$N], K=rep$value[data$N*2+1], K.sd=rep$sd[data$N*2+1], r=rep$value[data$N*2+2],r.sd=rep$sd[data$N*2+2],Q=rep$value[data$N*2+3],Q.sd=rep$sd[data$N*2+3],sigma=rep$value[data$N*2+4],sigma.sd=rep$sd[data$N*2+4],tau=rep$value[data$N*2+5],tau.sd=rep$sd[data$N*2+5],nll=rep$value[data$N*2+6],nll.sd=rep$sd[data$N*2+6])
      
      
    ###########################################################################
    # Results plots
    ###########################################################################
    plot(data$IOA,ylim=c(0,max(MLE$predI+MLE$predI.sd)),pch=16,col='red')
    with(MLE,lines(year,predI))
    with(MLE,lines(year,predI+predI.sd,lty=2))
    with(MLE,lines(year,predI-predI.sd,lty=2))

  }

}


