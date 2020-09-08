
fishery_model = function(  p, DS="stan", plotresults=TRUE, ... ) {

  if (0) {

    year.assessment=2016
    p = bio.snowcrab::load.environment( year.assessment=year.assessment)
    p$fishery_model = list()
    p$fishery_model$outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment )

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
        vector <lower=eps>[U] K;
        vector <lower=eps,upper=3>[U] r;
        vector <lower=eps,upper=3>[U] q;
        vector <lower=eps,upper=3>[U] qs;
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
        matrix[MN,U] pd;
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] P;
        matrix[MN,U] C;

        matrix[MN,U] F;
        matrix[M,U] TAC;


        // -------------------
        // annual production
         for(j in 1:U) {
           pd[1,j] = bm[2,j]- bm[1,j] + rem[1,j] ; // approximation
           for (i in 2:N ){
             pd[i,j] = (bm[i+1,j]- bm[i-1,j])/2 + rem[i,j] ; // linear interpolation cancels out the bm[i,j] term
           }
           for(i in N1:(MN-1)) {
             pd[i,j] = (bm[i+1,j]- bm[i-1,j])/2 + er * bm[i-1,j] ;  // linear interpolation cancels out the bm[i,j] term
           }
           pd[MN,j] = (bm[MN,j]- bm[(MN-1),j]) + er * bm[(MN-1),j]  ; // approximation
         }

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:3) {
           for (i in 1:N) {
             F[i,j] =  1.0 - rem[i,j] / bm[i,j]  ;
           }
         }
         // spring fishery
//         {
//           int j;
//           j=3;
//           F[1,j] =  1.0 - rem[1,j] / bm[1,j] ; // approximation
//           for (i in 2:N) {
//             F[i,j] =  1.0 - rem[i-1,j] / bm[i,j]  ;
 //          }
 //        }

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
      //    BX2MSY[j] = 1.0 - step( bm[N1,j]-0.25 ) ; // test if bm >= 1/2 bmY
      //    Bdrop[j]  = 1.0 - step( bm[N1,j]-bm[N,j] ) ; // test if bm(t) >= bm(t-1)
      //    Fcrash[j] = 4.0 * MSY[j] / exp(K[j]) ; // fishing mortality at which the stock will crash
        }

        // recaled estimates
         for(j in 1:U) {
           for(i in 1:MN) {
             B[i,j] = (bm[i,j] - rem[i,j]) * K[j] ;
             P[i,j] = pd[i,j]*K[j] ;
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



  if (DS=="stan_data" ) {
    sb = snowcrab_tsdata(p=p, assessment_years=2000:p$year.assessment)
    return(sb)
  }



  if (DS=="stan" ) {

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

}
