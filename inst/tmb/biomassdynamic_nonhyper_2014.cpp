#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
  
  // Data
  DATA_INTEGER(U);
  DATA_INTEGER(N);
  DATA_ARRAY(IOA);
  DATA_ARRAY(CAT);

 
  // Parameters:
  PARAMETER_VECTOR(log_sigma);       // process noise
  PARAMETER_VECTOR(log_tau);   // observation noise for commercial
  PARAMETER_VECTOR(q);                   // catchability for commercial size
  PARAMETER_VECTOR(log_r);                   // catchability for commercial size
  PARAMETER_VECTOR(log_K);                   // catchability for commercial size
  PARAMETER_ARRAY(log_P);        // unobserved state, commercial size
  
  vector<Type> sigma=exp(log_sigma);
  vector<Type> tau=exp(log_tau);
  vector<Type> r=exp(log_r);
  vector<Type> K=exp(log_K);
  matrix<Type> P=exp(log_P);

  // ############################################################################
  // Objective function: negative log likelihood
  Type nll = 0.0; // initialize neg loglik
  
  // -------------------  

  for (int j=1; j<U; j++) {
   
    for(int i=1; i<N; i++){
    
      // Dynamics equation (Pmed is the log of the predicted scaled biomass from the biomass dynamics)
      Pmed[i,j] = P[i-1,j] + r * P[i-1,j] * (1 - P[i-1,j]) - CAT[i-1,j]/K[j];
      
      // Process equation (Pmed is fit to the 'actual' scaled biomass (P) using the estimated process error (sigma) 
      nll -= dnorm(log(P[i,j]), log(Pmed[i,j]), sigma[j], true);
    }

    // Observation model: commercial size
    vector<Type> Omed(N);
    for(int i=0; i<NY; i++){

      // Observation equation (Omed is the log of the predicted index of abundance)
      Omed[i,j] = log(q[j] * K[j] * P[i,j]);      

      nll -= dnorm(log(IOA[i,j]), log(Omed[i,j]), tau[j], true);
    }
  }
    
  return(nll);
}




