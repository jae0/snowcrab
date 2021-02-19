
fishery_model = function(  p,  DS="logistic", assessment_years=2000:p$year.assessment, plotresults=TRUE, tag="default", areas=c("cfanorth", "cfasouth", "cfa4x"), ... ) {

#     sb = bio.snowcrab::fishery_model( DS="data_aggregated_timeseries", p=p, assessment_years=p$yrs )

  if (tag=="default") {
    if (exists("tag", p)) tag = p$tag 
  }

  
  if (0) {
    year.assessment=2016
    p = bio.snowcrab::load.environment( year.assessment=year.assessment)
    p$fishery_model = list()
    p$fishery_model$outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment )
  }


  if (DS=="logistic_parameters") {
    
    if (!exists("fishery_model", p)) p$fishery_model = list()
    
    if (!exists("method", p$fishery_model)) p$fishery_model$method = "stan"  # "jags", etc.
    if (!exists("outdir", p$fishery_model)) p$fishery_model$outdir = file.path( p$modeldir, "fishery_model_results"  )
    if (!exists("fnres", p$fishery_model)) p$fishery_model$fnres  = file.path( p$fishery_model$outdir, paste( "surplus.prod.mcmc", p$year.assessment, p$fishery_model$method, tag, "rdata", sep=".") )
    
    message( "Results will be saved to:", p$fishery_model$outdir)

    # observations
    if (!exists("standata", p$fishery_model)) p$fishery_model$standata = fishery_model( DS="data_aggregated_timeseries", p=p, assessment_years=p$yrs )     
      
    if (!exists("er", p$fishery_model$standata)) p$fishery_model$standata$er = 0.2  # target exploitation rate
    if (!exists("U", p$fishery_model$standata))  p$fishery_model$standata$U = ncol( p$fishery_model$standata$B)  # number of regions
    if (!exists("N", p$fishery_model$standata))  p$fishery_model$standata$N = nrow( p$fishery_model$standata$B)  # no years with data
    if (!exists("M", p$fishery_model$standata))  p$fishery_model$standata$M = 3 # no years for projections
    if (!exists("ty", p$fishery_model$standata)) p$fishery_model$standata$ty = which(p$assessment_years == 2004)  # index of the transition year (2004) between spring and fall surveys
    if (!exists("cfa4x", p$fishery_model$standata))  p$fishery_model$standata$cfa4x = 3 # column index of cfa4x
    if (!exists("eps",   p$fishery_model$standata))  p$fishery_model$standata$eps = 1e-6  # small non-zero number

    p$fishery_model$standata$missing = ifelse( is.finite(p$fishery_model$standata$IOA), 0, 1)
    p$fishery_model$standata$missing_n = colSums(p$fishery_model$standata$missing)
    p$fishery_model$standata$missing_ntot = sum(p$fishery_model$standata$missing_n)

    # this must be done last
    p$fishery_model$standata$IOA[ which(!is.finite(p$fishery_model$standata$IOA)) ] = 0 # reset NAs to 0 as stan does not take NAs
    p$fishery_model$standata$CAT[ which(!is.finite(p$fishery_model$standata$CAT)) ] = p$fishery_model$standata$eps  # remove NA's

    # priors
    if (!exists("Kmu", p$fishery_model$standata)) p$fishery_model$standata$Kmu =  c( 5, 50, 1 )
    if (!exists("rmu", p$fishery_model$standata)) p$fishery_model$standata$rmu = c(1, 1, 1)
    if (!exists("qmu", p$fishery_model$standata)) p$fishery_model$standata$qmu = c(1, 1, 1)
    if (!exists("Ksd", p$fishery_model$standata)) p$fishery_model$standata$Ksd =  c(0.5, 0.5, 0.5) * p$fishery_model$standata$Kmu  # c( 2, 20, 0.5)
    if (!exists("rsd", p$fishery_model$standata)) p$fishery_model$standata$rsd =  c(0.5, 0.5, 0.5) * p$fishery_model$standata$rmu  # rep( 0.3, 3)
    if (!exists("qsd", p$fishery_model$standata)) p$fishery_model$standata$qsd =  c(0.5, 0.5, 0.5) * p$fishery_model$standata$qmu  # rep( 0.3, 3)
    
    if (!exists("stancode", p$fishery_model )) p$fishery_model$stancode = fishery_model( p=p, DS="stan_surplus_production" )
    if (!exists("stancode_compiled", p$fishery_model )) p$fishery_model$stancode_compiled = rstan::stan_model( model_code=p$fishery_model$stancode )

    return(p)
  }


  if (DS=="stan_surplus_production") {
    return( "
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
        vector <lower=eps> [U] K;
        vector <lower=eps, upper=2> [U] r;
        vector <lower=eps, upper=2> [U] q;
        vector <lower=eps, upper=2> [U] qs;
        vector <lower=eps, upper=(1-eps)> [U] bosd;  // observation error
        vector <lower=eps, upper=(1-eps)> [U] bpsd;  // process error
        vector <lower=eps, upper=(1-eps)> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps> [M+N,U] bm;
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
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018

        for (j in 1:2) {
          Ymu[1,j]        = qs[j] * bm[1,j] - rem[1,j] ; // starting year approximation
          Ymu[2:(ty-1),j] = qs[j] * bm[2:(ty-1),j] - rem[1:(ty-2),j] ; //spring surveys
          Ymu[ty,j]       = q[j]  * bm[ty,j] - (rem[(ty-1),j] + rem[ty,j] )/2.0  ; // transition year .. approximation
          Ymu[(ty+1):N,j] = q[j]  * bm[(ty+1):N,j] - rem[(ty+1):N,j] ;   // fall surveys
        }
        {
          int k;
          k=3;
          Ymu[1,k]        = q[k] * bm[1,k]   - rem[1,k] ; // starting year approximation ymu[1991] = bm[1991]-rem[1991]
          Ymu[2:(ty-1),k] = q[k] * bm[2:(ty-1),k] - rem[2:(ty-1),k];
          Ymu[ty:N,k]     = q[k] * bm[ty:N,k] - rem[ty:N,k];
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
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;


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
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        matrix[M,U] TAC;


        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:3) {
           for (i in 1:N) {
             F[i,j] =  1.0 - rem[i,j] / bm[i,j]  ;
           }
         }
         for (j in 1:U) {
           for (i in N1:MN) {
             F[i,j] =  1.0 - er * bm[i-1,j] / bm[i,j]  ;
           }
           for (i in 1:MN) {
             F[i,j] =  -log( fmax( F[i,j], eps) )  ;
           }
         }

        // -------------------
        // parameter estimates for output

        for(j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
         for(j in 1:U) {
           for(i in 1:MN) {
             B[i,j] = (bm[i,j] - rem[i,j]) * K[j] ;
             C[i,j] = rem[i,j]*K[j] ;
           }
           for(i in 1:M) {
             TAC[i,j] = rem[N+i,j]*K[j] ;
           }
         }

      }
    "
    )
  }



  if (DS=="data_aggregated_timeseries" ) {


    cfanorth =  1 # column index
    cfasouth =  2 # column index
    cfa4x =  3 # column index

    landings = bio.snowcrab::snowcrab_landings_db()
      # NOTE:: message( "Fishing 'yr' for CFA 4X has been set to starting year:: 2001-2002 -> 2001, etc.")
      # year is year of capture
      # yr is "fishing year" relative to the assessment cycle
    landings = landings[ which (landings$cfa %in% c( "cfanorth", "cfasouth", "cfa4x" ) ) , ]
    L = tapply( landings$landings, INDEX=landings[,c("yr", "cfa")], FUN=sum, na.rm=T )
    nL = nrow(L)

    cfaall = tapply( landings$landings, INDEX=landings[,c("yr")], FUN=sum, na.rm=T )
    L = cbind( L, cfaall )
    L = L / 1000/1000  # convert to kt
    L[ !is.finite(L)] = 0

    L = as.data.frame( L[ match( assessment_years, rownames(L) ), areas ] )

    # biomass data: post-fishery biomass are determined by survey B)
    B = snowcrab.db(p=p, DS="carstm_output_timeseries"  )

    rownames(B) = B$yrs
    B = as.data.frame( B[ match( assessment_years, B$yrs ), areas ] )

    # cfa4x have had no estimates prior to 2004

    cfanorth.baddata = which( assessment_years <= 1997 )
    B[ cfanorth.baddata, cfanorth ] = NA

    cfasouth.baddata = which( assessment_years <= 1998 )
    B[ cfasouth.baddata, cfasouth ] = NA

    cfa.nodata =   which( assessment_years <= 2003 )
    B[ cfa.nodata , cfa4x ] = NA

    sb = list(
      IOA = as.matrix(B), # observed index of abundance
      CAT = as.matrix(L)  # catches  , assume 20% handling mortality and illegal landings
    )
    
    return(sb)
  }



  if (DS=="logistic_model" ) {

    library(rstan)
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())


    message( "Output location is: ", p$fishery_model$outdir )

    dir.create( p$fishery_model$outdir, recursive=T, showWarnings=F )

    f = rstan::sampling( p$fishery_model$stancode_compiled, data=p$fishery_model$standata, ... )
          # warmup = 200,          # number of warmup iterations per chain
          # control = list(adapt_delta = 0.9),
          # # refresh = 500,          # show progress every 'refresh' iterations
          # iter = 1000,            # total number of iterations per chain
          # chains = 5,             # number of Markov chains
          # cores = 5              # number of cores (using 2 just for the vignette)

    res = list( mcmc=rstan::extract(f), p=p)
    save(res, file=p$fishery_model$fnres, compress=TRUE)
    return(res)
  }


  if (DS=="logistic_samples" ) {
    res = NULL
    if (file.exists(p$fishery_model$fnres)) load(p$fishery_model$fnres)
    return(res)
  }


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



