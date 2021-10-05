// This model assumes two normal distributions (one for disease/one for not disease)
data {
   int N;                                // number of patients
   int Ntest;                            // number of tests
   vector[Ntest] y[N];                   // data
   int<lower=0> K_sites;                 // number of different sites (populations)
   int<lower=1,upper=K_sites> site[N];   // index of sites for each patient
   vector[2] prior_mu[Ntest];            // prior on mixture mean values
   vector[2] prior_sigma_mu[Ntest];      // prior on mixture mean values
   vector[2] beta_prior_prev[K_sites];
   int<lower=1,upper=2> cluster[Ntest,2];
}
parameters {
   real<lower=0,upper=1> prev[K_sites];   //mixing proportions
   ordered[2] mu[Ntest];
   vector<lower=0>[2] sigma[Ntest];
}

model {
   // **** Prior **** //
   for(t in 1:Ntest){
      for(i in 1:2){
         mu[t][i] ~ normal(prior_mu[t][i],prior_sigma_mu[t][i]);
         sigma[t][i] ~ exponential(2);
      }
   }
   for(ks in 1:K_sites){
      prev[ks] ~ beta(beta_prior_prev[ks][1],beta_prior_prev[ks][2]);
   }

   // **** Likelihood **** //
   for(n in 1:N){
      int k = site[n];
      real p1 = log(prev[k]);
      real p2 = log1m(prev[k]);
      for(t in 1:Ntest){
         p1 += normal_lpdf(y[n][t] | mu[t][cluster[t,1]], sigma[t][cluster[t,1]]);
         p2 += normal_lpdf(y[n][t] | mu[t][cluster[t,2]], sigma[t][cluster[t,2]]);
         //print("loop iteration: ",cluster[t,1], cluster[t,2])
      }
      target += log_sum_exp(p1,p2);
   }

}

generated quantities {
   real ps_SM[N];
   // This computes respective densities of each mixture component
   for(n in 1:N){
      int k = site[n];
      real p1 = prev[k];
      real p2 = 1-prev[k];
      for(t in 1:Ntest){
         p1 = p1*exp(normal_lpdf(y[n][t] | mu[t][cluster[t,1]], sigma[t][cluster[t,1]]));
         p2 = p2*exp(normal_lpdf(y[n][t] | mu[t][cluster[t,2]], sigma[t][cluster[t,2]]));
      }
      // normalise to get P_SM
      ps_SM[n] = p1/(p1+p2);
   }
}

