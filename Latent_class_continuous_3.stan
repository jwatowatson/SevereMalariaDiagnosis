//**** Three component model ****
// This model assumes a univariate normal distribution for *disease*
// and a mixture of two normals as the distribution for *not disease*
data {
   int N;                                // number of patients
   int Ntest;                            // number of tests
   vector[Ntest] y[N];                   // data
   int<lower=0> K_sites;                 // number of different sites (populations)
   int<lower=1,upper=K_sites> site[N];   // index of sites for each patient
   vector[3] prior_mu[Ntest];            // prior on mixture mean values
   vector[3] prior_sigma_mu[Ntest];      // prior on mixture mean values
   vector[3] alpha_prior[K_sites];       // prior on mixtures
   int<lower=1,upper=3> cluster[Ntest,3];
}

parameters {
   simplex[3] theta[K_sites];   //
   ordered[3] mu[Ntest];
   vector<lower=0>[3] sigma[Ntest];
}

model {
   // **** Prior **** //
   for(t in 1:Ntest){
      for(i in 1:3){
         mu[t][i] ~ normal(prior_mu[t][i],prior_sigma_mu[t][i]);
         sigma[t][i] ~ exponential(2);
      }
   }
   for(ks in 1:K_sites){
      theta[ks] ~ dirichlet(alpha_prior[ks]);
   }

   // **** Likelihood **** //
   for(n in 1:N){
      int k = site[n];
      vector[3] ps;
      // loop over each sub-distribution in the mixture
      for(i in 1:3){
         ps[i]=log(theta[k][i]);
         for(t in 1:Ntest){
            ps[i] += normal_lpdf(y[n][t] | mu[t][cluster[t,i]], sigma[t][cluster[t,i]]);
         }
      }
      target += log_sum_exp(ps);
   }

}

generated quantities {
   real<lower=0> ps_SM[N];
   vector<lower=0>[3] ps_comp[N];
   // This computes respective densities of each mixture component
   for(n in 1:N){
      int k = site[n];
      for(i in 1:3){
         ps_comp[n][i]=theta[k][i];
         for(t in 1:Ntest){
            ps_comp[n][i] *= exp(normal_lpdf(y[n][t] | mu[t][cluster[t,i]], sigma[t][cluster[t,i]]));
         }
      }
      // normalise to get P_SM
      ps_comp[n] = ps_comp[n]/sum(ps_comp[n]);
      ps_SM[n] = ps_comp[n][1];
   }
}
