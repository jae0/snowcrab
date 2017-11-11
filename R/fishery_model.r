
fishery_model = function(  p, DS="stan", plotresults=TRUE ) {

  if (0) {
    
    year.assessment=2016
    p = bio.snowcrab::load.environment( year.assessment=year.assessment)
    p$fishery_model = list()
    p$fishery_model$outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment )

  }

  message( "Output location is: ", p$fishery_model$outdir )

  dir.create( p$fishery_model$outdir, recursive=T, showWarnings=F )
 

  if (DS=="stan" ) {

    library(rstan)
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())

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
        matrix[MN,U] bmmu; // collator used to force positive values for lognormal
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


      }

      model {

        // -------------------  
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        q ~ normal( qmu, qsd )  ;
        qs ~ normal( qmu, qsd )  ;
        b0 ~ beta( 8, 2 ) ; // starting b prior to first catch event 
        bosd ~ cauchy( 0, 0.5 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.5 ) ;


        // -------------------  
        // biomass observation model 
        for (j in 1:U) {  
          log(Y[1:N,j]) ~ normal( log(Ymu[1:N,j]), bosd[j] ) ;  // stan thinks Y is being transformed due to attempt to impute missing values .. ignore
        }
 

        // -------------------  
        // biomass process model 
        for (j in 1:U) {
          log(bm[1:MN,j]) ~ normal( log(bmmu[1:MN,j]), bpsd[j] ) ;     
        }

        // could have used lognormal but this parameterization is 10X faster and more stable
        target += - log(fabs(Y));  // required due to log transf above
        target += - log(fabs(bm));

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

    f = sampling(stanmodel, data=sb, chains=4, iter=10000, warmup = 2000,
      control = list(adapt_delta = 0.9, max_treedepth=15) )
          # warmup = 200,          # number of warmup iterations per chain
          # control = list(adapt_delta = 0.9),
          # # refresh = 500,          # show progress every 'refresh' iterations
          # iter = 1000,            # total number of iterations per chain
          # chains = 5,             # number of Markov chains
          # cores = 5              # number of cores (using 2 just for the vignette)

   
    res = list( mcmc=extract(f), sb=sb, p=p)
    save(res, file=p$fishery_model$fnres, compress=T)
    return(res)

      if (0) {

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
 
  }



  # --------------------------
  


  if (DS=="jags") {

    warning( "Jags method is deprecated. Use stan.")

    require(rjags)
    rjags::load.module("dic")
    rjags::load.module("glm")

    if (!exists("fishery_model", p)) p$fishery_model = list()
    if (!exists("n.adapt", p$fishery_model)) p$fishery_model$n.adapt  = 5000  # burn-in  .. 4000 is enough for the full model 
    if (!exists("n.iter ", p$fishery_model)) p$fishery_model$n.iter   = 10000  #  n.iter = 30000,
    if (!exists("n.chains", p$fishery_model)) p$fishery_model$n.chains  = 8  #  n.chains = 8 ,
    if (!exists("n.thin", p$fishery_model)) p$fishery_model$n.thin  = 100  # high autocorrelations 
    
    n.iter.total = p$fishery_model$n.iter * p$fishery_model$n.thin

    sb = biomass.summary.db(p=p, DS="surplusproduction" )

    sb$tomonitor = c( "r", "K", "q", "qs", "r.mu", "r.sd", "b","bp.sd", "bo.sd","b0", "b0.sd", "rem", "rem.sd", "rem.mu","REM", "MSY", "BMSY", "FMSY", "Fcrash", "Bdrop", "BX2MSY", "F", "TAC",  "C", "P", "B" )

    sb$jagsmodelname = "biomassdynamic_nonhyper_2016.bugs"
    sb$jagsmodel = 
    "
model {

  # -------------------  
  # define some marginally informative variance priors using CV's (coefficients of variation) as a simple approach
  # uniform distribution seems most stable .. too uninformative?
  # NOTE: eps = a small number non-zero number (essentially equivalent to zero but used to prevent infinity values)
  # uninformative CV's associated with process (bp.) and observation (bo.) errors 
 
  # NOTE for lognormals: CV = sqrt(exp( SD ^2) - 1)  and CV ~ SD where SD ~ < 0.5
  # of SD = sqrt(log(1+CV^2)) and therefore, in terms of precision:
  # TAU = 1/SD^2 = 1/log(1+CV^2)


  for (j in 1:U) {
    r[j] ~ dnorm( r.mu[j], pow( r.sd[j], -2 ) ) 
  }


  for (j in 1:U) {
    #q.mu[j]  ~ dunif( q.min, q.max ) 
    #q.sd[j]  ~ dunif( q.mu[j] * cv.normal.min, q.mu[j] *cv.normal.max )  # catchability coefficient (normal scale)
    q[j] ~ dnorm( q.mu[j], pow( q.sd[j], -2 ) )  T(q.min, q.max)
    #q[j] ~ dunif( q.min[j] , q.max[j] )
  }


  for (j in 1:U) {
    K[j] ~ dlnorm( K.mu[j], pow( K.sd[j], -1 )) # T(K.min[j], K.max[j] )
  }


  # -------------------  
  # removals (catch) observation model, standardized to K (assuming no errors in observation of catch!)
    for (j in 1:U) {
      for (i in 1:N){
        rem[i,j] <- CAT[i,j]/K[j]
      }
    }


  # -------------------  
  # biomass observation model 
  #   This is slightly complicated because a fall / spring survey correction is required:
  #   B represents the total fishable biomass available in fishing year y
  #     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t) 
  #     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1) 
  #   this is conceptualized in the following time line: 
  #     '|' == start/end of each new fishing year
  #     Sf = Survey in fall
  #     Ss = Survey in spring
  #     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...

    for (j in 1:(U)) {
      #bo.tau[j] ~ dunif( pow( log( 1 + pow( cv.lognormal.max, 2) ), -1 ), pow( log( 1 + pow( cv.lognormal.min, 2) ), -1 ) )  # min/max inverted because it is an inverse scale
      #bo.sd[j] ~ dunif( bo.min[j], bo.max[j] )  
      bo.sd[j] ~ dlnorm(bo.mup[j],pow(bo.sdp[j],-1))
      bo.tau[j]  <- pow( bo.sd[j], -2 )
    }

    for (j in 1:(U-1)) {
      # spring surveys from 1998 to 2003
      IOA[1,j] ~ dlnorm( log( max( q[j] * K[j] * (bm[1,j] - rem[1,j]) , eps)), bo.tau[j] )  # approximation
      for (i in 2:(ty-1)) { 
        IOA[i,j] ~ dlnorm( log( max( q[j] * K[j] * (bm[i,j]- rem[(i-1),j]), eps)), bo.tau[j] )  ;
      }
      # transition year
      IOA[ty,j] ~ dlnorm( log( max( q[j] * K[j] * (bm[ty,j] - (rem[(ty-1),j] + rem[ty,j] )/2 ), eps)), bo.tau[j] ) ;  # approximation
      # fall surveys    
      for (i in (ty+1):N) {
        IOA[i,j] ~ dlnorm( log( max( q[j] * K[j] * (bm[i,j] - rem[i,j]), eps)), bo.tau[j] ) ;
      }
    }

    # Cfa 4X -- fall/winter fishery
    # assume similar to a spring fishery but no need for separate q's
    #    Btot(t) = Bsurveyed(t)+ removals(t-1)
    #    NOTE: year designation in 4X is for the terminal year: ie. 2001-2002 => 2002
    
    IOA[1,cfa4x] ~ dlnorm( log( max( q[cfa4x] * K[cfa4x] * (bm[1,cfa4x] - rem[1,cfa4x]), eps)), bo.tau[cfa4x] ) ;  # approximation
    for (i in 2:N) { 
      IOA[i,cfa4x] ~ dlnorm( log( max( q[cfa4x] * K[cfa4x] * (bm[i,cfa4x]- rem[(i-1),cfa4x]), eps)), bo.tau[cfa4x] ) ;
    }



  # -------------------  
  # biomass process model 
        
    for (j in 1:U) {
      #bp.tau[j] ~ dunif( pow( log( 1 + pow( cv.lognormal.max, 2) ), -1 ), pow( log( 1 + pow( cv.lognormal.min, 2) ), -1 ) )  
      bp.sd[j] ~ dunif( bp.min[j], bp.max[j] )  
      #bp.sd[j]  <- pow( sqrt(bp.tau[j]), -1 )
      # bp.sd[j] ~ dlnorm(bp.mup[j],pow(bp.sdp[j],-1))
      bp.tau[j]  <- pow( bp.sd[j], -2 )
    }

    for(j in 1:U) {
      b0[j] ~ dunif( b0.min[j], b0.max[j] ) # starting b prior to first catch event 
      bm[1,j] ~ dlnorm( log( max( b0[j], eps)), bp.tau[j] ) T(b.min, b.max ) ;  # biomass at first year   
      for(i in 2:(N+M)) {
        bm[i,j] ~ dlnorm( log( max(bm[i-1,j]*( 1 + r[j]*(1-bm[i-1,j])) - rem[i-1,j] , eps)), bp.tau[j] ) T(b.min, b.max) ;
      }
      
      # forecasts
      for(i in 1:M) {
        rem[N+i,j] <- er*bm[N+i-1,j]
      }
    }


  # -------------------  
  # monitoring nodes and parameter estimates for output
    for(j in 1:U) {
      Bdrop[j]  <- 1 - step( bm[N+1,j]-bm[N,j] ) ; # test if bm(t) >= bm(t-1) 
      BX2MSY[j] <- 1 - step( bm[N+1,j]-0.25 ) ; # test if bm >= 1/2 bmY
      MSY[j]    <- r[j]* exp(K[j]) / 4  # maximum height of of the latent productivity (yield)
      BMSY[j]   <- exp(K[j])/2  # biomass at MSY
      FMSY[j]   <- 2 * MSY[j] / exp(K[j]) # fishing mortality at MSY
      Fcrash[j] <- 4 * MSY[j] / exp(K[j]) # fishing mortality at which the stock will crash
    }


    # -------------------  
    # fishing mortality
    # force first year estimate assuming catches in year 0 to be similar to year 1 
    for(j in 1:U) {
      for(i in 1:N) {
        F[i,j] <- -log( max(1 - rem[i,j] / bm[i,j], eps))  
      }
      for(i in (N+1):(N+M)) {
        F[i,j] <- -log( max(1 - er * bm[i-1,j] / bm[i,j], eps)) 
      }
    }


    # -------------------  
    # annual production
    for(j in 1:U) {
      pd[1,j] <- bm[2,j]- bm[1,j] + rem[1,j] # approximation
      for (i in 2:(N) ){
        pd[i,j] <- (bm[i+1,j]- bm[i-1,j])/2 + rem[i,j]  # linear interpolation cancels out the bm[i,j] term
      }
      for(i in (N+1):(N+M-1)) {
        pd[i,j] <- (bm[i+1,j]- bm[i-1,j])/2 + er * bm[i-1,j]   # linear interpolation cancels out the bm[i,j] term
      }
      pd[(N+M),j] <- (bm[(N+M),j]- bm[(N+M-1),j]) + er * bm[(N+M-1),j]   # approximation
    }
  


    # -------------------  
    # recaled estimates
  
    for(j in 1:U) {
      for(i in 1:(N+M)) {
        B[i,j] <- (bm[i,j] - rem[i,j]) * K[j]
        P[i,j] <- pd[i,j]*K[j]
        C[i,j] <- rem[i,j]*K[j]
      }
      for(i in 1:M) {
        TAC[i,j] <- rem[N+i,j]*K[j]
      }
    }

}

    "

    fn = tempfile()
    cat( sb$jagsmodel, file=fn )

    m = jags.model( file=fn, data=sb, n.chains=p$fishery_model$n.chains, n.adapt=p$fishery_model$n.adapt ) # recruitment + spring/summer q's + all observed CVs
    
    tomonitor = intersect( variable.names (m), sb$tomonitor )
    
    y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.total, thin=p$fishery_model$n.thin) # sample from posterior
    
    jags2stan = function(x){
      # merge chains together
      print(attr(x, "varname"))
      xdim = dim(x)
      chain = length(xdim)
      iter = chain-1
      oorder = c(iter, chain, 1:(iter-1))
      o = aperm( x,  oorder)
      odim = dim(o)
      dim(o) =c(odim[1]*odim[2], odim[3:length(odim)])
      return( o)
    }
    out = lapply( y, jags2stan )

    res = list( mcmc=out, sb=sb, p=p)
    save(res, file=p$fishery_model$fnres, compress=T)

    return(res)


      if (0) {

        dic.samples(m, n.iter=p$fishery_model$n.iter ) # pDIC

        graphics.off() ; plot.new(); layout( matrix(c(1,2,3), 3, 1 )); par(mar = c(5, 4, 0, 2))
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
        #  update(m, n.iter=p$fishery_model$n.iter ) # above seems enough for convergence but a few more to be sure

        # determine autocorrelation thinning
        # y = coda.samples(m, variable.names=c("K", "r", "q"), n.iter=20000, thin=10) # sample from posterior
        # autocorr.plot(y)   # about 10 to 20 required
        # plot(y, ask=T)
        # autocorr(y, lags = c(0, 1, 5, 10, 50), relative=TRUE)

        # final sampling from the posteriors
        #  y = jags.samples(m, variable.names=tomonitor, n.iter=10000, thin=20) # sample from posterior
        y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.total, thin=p$fishery_model$n.thin) # sample from posterior
        
        

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


        plot.new()
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


  # --------------------------


  if (DS=="LaplacesDemon") {

    warning( "LaplacesDemon method is not yet complete")

    require(LaplacesDemonCpp)
      

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

      if (is.null( Data$Missing )) {
        mO = which( !is.finite(Data$O))
        mr = which( !is.finite(Data$removals) )
        Data$Missing = list(
          O = mO,
          nO = length(mO) ,
          removals = mr , 
          nremovals = length( mr ),
          n = length( mO) + length(mr) 
        )
      }

      if (is.null( Data$PGF )) {
        # Parameter Generating Function
        Data$PGF = function(Data) {
          r=rnorm(Data$nregions, Data$r0, sd=Data$cv )
          K=rnorm(1, log(Data$K0), sd=Data$cv ) # log scale
          q=rnorm(1, Data$q0, sd=Data$cv )
          S_sd=runif( 1 )
          O_sd=runif( 1 )
          S=rnorm( Data$N, log(Data$S0), sd=Data$cv ) # log scale 
          S0=rnorm(1, log(Data$S0), sd=Data$cv ) # log scale
          out = c(r, K, q, S_sd, O_sd, S, S0)
          if (Data$Missing$nO > 0) {
            Omissing = rnorm( Data$Missing$nO, log(Data$Omissing0), sd=Data$cv )
            out = c( out, Omissing )
          }
          if (Data$Missing$nremovals > 0) {
            removalsmissing = rnorm( Data$Missing$nremovals, log(Data$removalsmissing0), sd=Data$cv )
            out = c( out, removalsmissing )
          }
          return( out  )
        }
        Data$PGF = compiler::cmpfun(Data$PGF)
      }

      if (is.null( Data$parm.names )) {
        # paramater names and initial values
        pnames = list(
          r = Data$r0,
          K = Data$K0,
          q = Data$q0,
          S_sd = Data$cv,
          O_sd = Data$cv,
          S = rep( Data$S0, Data$N),
          S0 = Data$S0
        )
        if (Data$Missing$nO > 0) pnames$Omissing = rep( Data$Omissing0, Data$Missing$nO )
        if (Data$Missing$nremovals > 0)  pnames$removalsmissing = rep( Data$removalsmissing0, Data$Missing$nremovals )
        Data$parm.names = as.parm.names( pnames )
      }

      if (is.null( Data$idx )) {
        # index position of paramaters
        Data$idx = list(
          r=grep("\\<r\\>", Data$parm.names),
          K=grep("\\<K\\>", Data$parm.names),
          q=grep("\\<q\\>", Data$parm.names),
          S_sd=grep("\\<S_sd\\>", Data$parm.names),
          O_sd=grep("\\<O_sd\\>", Data$parm.names),
          S=grep("\\<S\\>", Data$parm.names),
          S0=grep("\\<S0\\>", Data$parm.names)
        )
        if (Data$Missing$nO > 0) Data$idx$Omissing = grep("\\<Omissing\\>", Data$parm.names)
        if (Data$Missing$nremovals > 0)  Data$idx$removalsmissing = grep("\\<removalsmissing\\>", Data$parm.names)
      }


      if (is.null( Data$mon.names )) {
        Data$mon.names = c("LP", "r", "K", "q" )
        # Data$mon.names = c("LP", "r", "K", "q", paste0("S",1:Data$N), paste0("AR",1:(Data$N-1) ) )
      }


      if (is.null( Data$Model )) {
        Data$Model = function(parm, Data) {
          
          Spred = Opred = rep(0, Data$N )   # initialize a few storage vectors

          # constraints:
          parm[Data$idx$q] = LaplacesDemonCpp::interval( parm[Data$idx$q], Data$eps, 2 );
          parm[Data$idx$r] = LaplacesDemonCpp::interval( parm[Data$idx$r], Data$eps, 2 );
          parm[Data$idx$K] = LaplacesDemonCpp::interval( parm[Data$idx$K], log(Data$eps), log(Data$K0*2) );
          parm[Data$idx$S0] = LaplacesDemonCpp::interval( parm[Data$idx$S0], log(Data$eps), log(Data$smax) );
          parm[Data$idx$S] = LaplacesDemonCpp::interval( parm[Data$idx$S], log(Data$eps), log(Data$smax) );

          parm[Data$idx$O_sd] = LaplacesDemonCpp::interval( parm[Data$idx$O_sd], Data$eps, 1)
          parm[Data$idx$S_sd] = LaplacesDemonCpp::interval( parm[Data$idx$S_sd], Data$eps, 1)

          # these are SD's on log scale (0,1) is sufficient
          # continue with SD priors as these are now on correct scale
          
          loglik = c()
          loglik[Data$parm.names] = 0
          loglik[Data$idx$q] = LaplacesDemonCpp::dhalfcauchy( parm[Data$idx$q], Data$cv,  TRUE ) ;
          loglik[Data$idx$r] = dnorm( parm[Data$idx$r], Data$r0, Data$cv, TRUE ) ;
          loglik[Data$idx$K] = dnorm( parm[Data$idx$K], log(Data$K0), Data$cv, TRUE ) ;
          loglik[Data$idx$S0] = dnorm( parm[Data$idx$S0], log(Data$S0), Data$cv, TRUE ) ;
          loglik[Data$idx$O_sd] = LaplacesDemonCpp::dhalfcauchy( parm[Data$idx$O_sd], Data$cv, TRUE );
          loglik[Data$idx$S_sd] = LaplacesDemonCpp::dhalfcauchy( parm[Data$idx$S_sd], Data$cv, TRUE );

          q = parm[Data$idx$q]
          r = parm[Data$idx$r]
          K = exp( parm[Data$idx$K] )
          Spred[1] = exp( parm[Data$idx$S0] );
          S = exp( parm[Data$idx$S] ) ;
          O_sd = parm[Data$idx$O_sd]
          S_sd = parm[Data$idx$S_sd]

          R = Data$removals/K ;
          if (exists( "Missing", Data ) ) {
            if (Data$Missing$nremovals > 0) {
              R[Data$Missing$removals ] = rlnorm( Data$Missing$nremovals, 
                parm[Data$idx$removalsmissing], sd=Data$cv )
              R = .Internal(pmax(na.rm=TRUE, R, Data$eps )) ;
              loglik[Data$idx$removalsmissing] = dnorm( parm[Data$idx$removalsmissing], 
                log(R[Data$Missing$removals ]), Data$cv )
            }
          }

          AR = 1.0 + r*(1-S[Data$jj]) ;  # autocorelation component
          Spred[ Data$ii ] = S[Data$jj] * AR - R[Data$jj] ;  # simple logistic
          Spred = .Internal(pmax(na.rm=TRUE, Spred, Data$eps )) ;
          loglik[Data$idx$S] = dnorm( parm[Data$idx$S], log(Spred), S_sd, TRUE );

          # Likelihoods for observation model
          i = 1 ;                   Opred[i] = S[i] - R[i] ;
          i = 2:(Data$ty-1);        Opred[i] = S[i] - R[i-1] ;
          i = Data$ty;              Opred[i] = S[i] - (R[i-1] + R[i])/2 ;
          i = (Data$ty+1):Data$N ;  Opred[i] = S[i] - R[i] ;
          Opred = .Internal(pmax( na.rm=TRUE, K*q*Opred, Data$eps ) ) ;
          if (exists( "Missing", Data ) ) {
            if (Data$Missing$nO > 0) {
              Data$O[Data$Missing$O ] = rlnorm( Data$Missing$nO, parm[Data$idx$Omissing], sd=Data$cv )
              Data$O[Data$Missing$O ] = .Internal(pmax(na.rm=TRUE, Data$O[Data$Missing$O ], Data$eps )) ;
              loglik[Data$idx$Omissing] = dnorm( parm[Data$idx$Omissing], log(Data$O[Data$Missing$O ]), Data$cv )
            }
          }
          ll_obs = dnorm( log(Data$O), log(Opred), O_sd, TRUE )

          Smon = Rmon = ER = F = B = C = rep(0, Data$MN )
          
          Smon[1:Data$N] = S
          Rmon[1:Data$N] = R
          for( i in (Data$N+1):Data$MN ){
            Smon[i] = Smon[i-1] * (1.0 + r*(1.0-Smon[i-1]) ) - Rmon[i-1]
            Rmon[i] = Smon[i-1] * Data$er 
          }
          Smon = .Internal(pmax( na.rm=TRUE, Smon, Data$eps ) )
          Rmon = .Internal(pmax( na.rm=TRUE, Rmon, Data$eps ) )
          
          # monitoring
          ER = Rmon / Smon ;
          B = Smon*K
          C = Rmon*K
          F = -log( 1 - ER) ; # fishing mortality
                
          LL = sum( ll_obs )  # log likelihood (of the data)
          LP = LL + sum( loglik )  # log posterior
          out = list( LP=LP, Dev=-2*LL, Monitor=c(LP, r, K, q), yhat=S*K, parm=parm )
          return( out  )
        }
        Data$Model.ML  = compiler::cmpfun( function(...) (Data$Model(...)$Dev / 2) )  # i.e. - log likelihood
        Data$Model.PML = compiler::cmpfun( function(...) (- Data$Model(...)$LP) ) #i.e., - log posterior 
        Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed
      }

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

    warning( "TMB method is not yet complete")

    require(stmenv)
      

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

  # -----------------------

  if (DS=="delay.difference_jags") {

    warning( "This is not yet complete .. will eventually convert to Stan for speed and stability ..")    

    require(rjags)
    rjags::load.module("dic")
    rjags::load.module("glm")

    bioLibrary( "bio.models")

    ###  all data follow this sequence: c("cfanorth", "cfasouth", "cfa4x")
    redo.data=F
    if (redo.data) {
      biomass.summary.db("complete.redo", p=p )
    }
    res = biomass.summary.db( p=p )


    sb = list(
      FB0x = c(0.8, 0.6, 0.1),  # priors of the starting fishable biomass : N,S,4X
      q0x = c(1, 1, 1),  # priors of "catchability"
      K0x = c(4, 50, 1 ),  # mean carrying capacity estimate
      IOA = as.matrix(res$B), # observed index of fishable biomass
      CAT = as.matrix(res$L) , # observed landings
      REC = as.matrix(res$R) , # observed recruitment
      qREC0 = 1, # q correction for recruitment
      sigma = 0.3, # default upper limit in error (CV) of data
      sigmaB = 0.3, # default upper limit in error (CV) of FB data
      cv = 0.5, # ...
      er = 0.2,  # target exploitation rate
      MFB = 0.8, # prior estimate of natural mortality for FB
      MFB.sd = 0.8, # prior estimate of SD in natural mortality for FB
      MREC = 0.8, # prior estimate of natural mortality for FB
      MREC.sd = 0.8, # prior estimate of SD in natural mortality for FB
      brodycoef = 2,  # proportional increase in weight with each age
      brodysd0 = 0.2,  # proportional increase in weight with each age
      N = nrow( res$B) , # no years with data
      M = 5, # no years for projections
      R = ncol( res$B),  # number of regions
      ty=7,  # index of the transition year (2004) between spring and fall surveys
      cfa4x=3, # index of cfa4x
      eps = 1e-4  # small non-zero number
    )

    # MCMC/Gibbs parameters
    n.adapt = 5000 # burn-in
    n.chains = 3
    n.thin = 100
    n.iter = 10000


    if (modelrun == "simple" ) {
      ## model 1 --- simple delay difference with observation and process error
      m = jags.model( file=fishery.model.jags ( DS="delay.difference" ), data=sb, n.chains=n.chains, n.adapt=n.adapt )
      coef(m)
      tomonitor = c("FB","K", "REC", "Znat", "Ztot", "q.biomass", "qREC", "Catch", "Cr" )
      dic.samples(m, n.iter=n.iter ) # pDIC
      fnres = file.path( project.datadirectory("bio.snowcrab"), "output", "delaydifference.mcmc.simple.rdata" )

    }


    if (modelrun == "simple.illegal" ) {
      ## model 2 --- simple delay difference with observation and process error and illegal landings
      m = jags.model( file=fishery.model.jags ( DS="delay.difference.illegal" ), data=sb, n.chains=n.chains, n.adapt=n.adapt )
      coef(m)
      tomonitor =  c("FB","K", "REC", "Znat", "Ztot", "q.biomass", "qREC", "Catch", "Cr" )
      dic.samples(m, n.iter=n.iter ) # pDIC
      fnres = file.path( project.datadirectory("bio.snowcrab"), "output", "delaydifference.mcmc.simple.illegal.rdata" )

    }


    tomonitor = intersect( variable.names (m), tomonitor )


    # convergence testing -- by 1000 to 1500 convergence observed by Gelman shrink factor diagnostic
    convergence.test = F
    if (convergence.test) {
      y = jags.samples(m, variable.names=tomonitor, n.iter=10000, thin=10)
      gelman.plot(y[["Ztot"]])
      gelman.plot(y[["K"]])
      gelman.plot(y[["q.biomass"]])
      gelman.plot(y[["qREC"]])
      gelman.plot(y[["sd.K"]])
      gelman.plot(y[["sd.o"]])
      gelman.plot(y[["sd.p"]])
      geweke.plot(y[["r"]])
    }


    # autocorrelation thinning
      y = coda.samples(m, variable.names=c("K", "r", "q"), n.iter=10000, thin=10) # sample from posterior
      autocorr.plot(y)
      # plot(y, ask=T)
      # autocorr(y, lags = c(0, 1, 5, 10, 50), relative=TRUE)


    # update if not yet converged
      update(m, n.iter=n.iter ) # above seems enough for convergence but a few more to be sure


    # final sampling from the posteriors
      n.iter.final = n.iter * n.thin
      # n.iter.final = n.iter
      y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.final, thin=n.thin) # sample from posterior


    # save(y, file=fnres, compress=T)
    # load( fnres )


    # Figures
      dir.output = file.path( dirname(p$ofname), "figures", "bugs")
      dir.create( dir.output, recursive=T, showWarnings=F )


      # frequency density of key parameters
      figure.bugs( "K", y=y, fn=file.path(dir.output, "K.density.png" ) )
      figure.bugs( "r", y=y, fn=file.path(dir.output, "r.density.png" ) )
      figure.bugs( "q", y=y, fn=file.path(dir.output, "q.density.png" ) )
      figure.bugs( "BMSY", y=y, fn=file.path(dir.output, "BMSY.density.png" ) )

      # timeseries
      figure.bugs( type="timeseries", var="biomass", y=y, fn=file.path(dir.output, "biomass.timeseries.png" ) )
      figure.bugs( type="timeseries", var="fishingmortality", y=y, fn=file.path(dir.output, "fishingmortality.timeseries.png" ) )

      # Harvest control rules
      figure.bugs( type="hcr", var="default", y=y, fn=file.path(dir.output, "hcr.default.png" ) )
      figure.bugs( type="hcr", var="simple", y=y, fn=file.path(dir.output, "hcr.simple.png" ) )

      # diagnostics
      figure.bugs( type="diagnostics.production", y=y, fn=file.path(dir.output, "diagnostic.production.png" ) )
      figure.bugs( type="diagnostics....", y=y, fn=file.path(dir.output, "diagnostic. ... .png" ) )




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

        # densities of F in previous year
        for (i in 1:3) plot(density( y$F[ndata-1,i,,] ), main="")
        qs = apply( y$F[ndata-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
        qs

        # F for table ---
        summary(y$F, median)




    # B timeseries with error
    graphics.off()
    layout( matrix(c(1,2,3), 3, 1 ))
    par(mar = c(5, 4, 0, 2))

    vv = "B"
    Xm = jags.extract( y, vv, mean )
    X = y[[vv]]
    nR = dim(X)[2]
    for (r in 1:nR) {
      nn = create.histograms.yearly( X[,r,,] )
      nn = nn[ 1:(nrow(nn)-1) ,]  # truncate larges and smallest data points as they are not linear increments
      nnr = as.numeric(rownames(nn) )
      if (r==2) nn = nn[ which( nnr < 100 ) ,]
      if (r==3) nn = nn[ which( nnr < 2 ) ,]
      colnames(nn) = yrs
      plot.histograms.temporal( nn, overlaydata=Xm[,r], xlab="Year", ylab="Biomass; kt", barscale=.9, barcol="gray" )
    }



    # F timeseries with error
    graphics.off()
    layout( matrix(c(1,2,3), 3, 1 ))
    par(mar = c(5, 4, 0, 2))

    vv = "F"
    Xm = jags.extract( y, vv, mean )
    X = y[[vv]]
    X[ X>1] = 1
    X = X[ 1:length(yrs0), ,, ]

    nR = dim(X)[2]
    for (r in 1:nR) {
      nn = create.histograms.yearly( X[,r,,] )
      nn = nn[ 2:(nrow(nn)-1) ,]  # truncate larges and smallest data points as they are not linear increments
      uu = scale( nn, center=FALSE, scale=colSums(nn) )
      colnames(nn) = yrs0
      plot.histograms.temporal( nn, overlaydata=Xm[1:length(yrs0),r], xlab="Year", ylab="Fishing mortality", barscale=0.9, barcol="gray" )
    }



  }




}


