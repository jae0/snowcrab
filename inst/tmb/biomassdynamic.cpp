
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
 
  // Data
  DATA_INTEGER(N);
  DATA_VECTOR(IOA);
  DATA_VECTOR(CAT);
 
  // Parameters:
  PARAMETER(log_sigmap);       // process noise
  PARAMETER(log_sigmao);   // observation noise for commercial
  PARAMETER(log_Q);                   // catchability for commercial size
  PARAMETER(log_r);                   // catchability for commercial size
  PARAMETER(log_K);                   // catchability for commercial size
  PARAMETER_VECTOR(log_P);        // unobserved state, commercial size
  
  Type sigmap=exp(log_sigmap);
  Type sigmao=exp(log_sigmao);
  Type Q=exp(log_Q);
  Type r=exp(log_r);
  Type K=exp(log_K);
  vector<Type> P=exp(log_P);

  vector<Type> Ppred(N);
  vector<Type> Ipred(N);
  vector<Type> B(N);

  // Objective function: negative log likelihood
  Type nll = 0.0; // initialize neg loglik
  
  // ############################################################################

   
  for(int i=1; i<N; i++){
  
    // Dynamics equation (Ppred is the log of the predicted scaled biomass from the biomass dynamics)
    Ppred(i) = P(i-1) + r * P(i-1) * (1 - P(i-1)) - CAT(i-1)/K;
    
    // Process equation (Ppred is fit to the 'actual' scaled biomass (P) using the estimated process error (sigma) 
    nll -= dnorm(P(i), Ppred(i), Ppred(i)*sigmap, true);
  }

  // Observation model: commercial size
  for(int i=0; i<N; i++){

    // Observation equation (Omed is the log of the predicted index of abundance)
    B(i) = P(i) * K;
    Ipred(i) = Q * B(i);      

    nll -= dnorm(IOA(i), Ipred(i), Ipred(i)*sigmao, true);
  }
  

  // ############################################################################
  

  // Report
  REPORT(B);
  REPORT(Ipred);
  REPORT(K);
  REPORT(r);
  REPORT(Q);
  REPORT(sigmap);
  REPORT(sigmao);
  REPORT(nll);
  
  ADREPORT(B);
  ADREPORT(Ipred);
  ADREPORT(K);
  ADREPORT(r);
  ADREPORT(Q);
  ADREPORT(sigmap);
  ADREPORT(sigmao);
  ADREPORT(nll);

    
  return nll;
}




