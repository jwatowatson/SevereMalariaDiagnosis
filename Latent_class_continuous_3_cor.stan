// This model assumes three bivariate distributions (one for disease/one for not disease)
// Component 1: not SM
// Component 2: not SM
// Component 3: SM

data {
   int N;                                // number of patients
   int Ntest;                            // number of tests
   vector[Ntest] y[N];                   // data
   vector[3] prior_mu[Ntest];
   vector<lower=0>[3] prior_mu_sigma[Ntest];
   vector<lower=0>[Ntest] prior_sigma_vals[3];
   int<lower=0> K_sites;                 // number of different sites (populations)
   int<lower=1,upper=K_sites> site[N];   // index of sites for each patient
   vector[3] dirichlet_prior_prev[K_sites];
   real<lower=0> nu;
}

parameters {
   simplex[3] prev_comp[K_sites];   //mixing proportions
   ordered[3] mu_test[Ntest];
   vector<lower=0>[Ntest] sigma[3];
   cholesky_factor_corr[Ntest] Lcorr[3];// cholesky factor (L_u matrix for R)
}

transformed parameters{
   corr_matrix[Ntest] R[3]; // correlation matrix
   cov_matrix[Ntest] Sigma[3]; // VCV matrix
   vector[Ntest] mu[3];

   for(t in 1:3){
      for(i in 1:Ntest){
         mu[t][i] = mu_test[i][t];
      }
   }

   for(i in 1:3){
      R[i] = multiply_lower_tri_self_transpose(Lcorr[i]);
      Sigma[i] = quad_form_diag(R[i], sigma[i]);
   }
}

model {
   // **** Prior **** //
   for(t in 1:3){
      for(i in 1:Ntest){
         sigma[t][i] ~ normal(prior_sigma_vals[t][i], .5);
      }
      Lcorr[t] ~ lkj_corr_cholesky(nu); // prior for cholesky factor
   }
   for(ks in 1:K_sites){
      prev_comp[ks] ~ dirichlet(dirichlet_prior_prev[ks]);
   }
   for(i in 1:Ntest){
      mu_test[i] ~ normal(prior_mu[i], prior_mu_sigma[i]);
   }
   // **** Likelihood **** //
   for(n in 1:N){

      int k = site[n];
      vector[3] ps;
      // groups 1 & 2 are not SM; group 3 is SM
      ps[1] = log(prev_comp[k][1])+multi_normal_lpdf(y[n] | mu[1], Sigma[1]);
      ps[2] = log(prev_comp[k][2])+multi_normal_lpdf(y[n] | mu[2], Sigma[2]);
      ps[3] = log(prev_comp[k][3])+multi_normal_lpdf(y[n] | mu[3], Sigma[3]);

      target += log_sum_exp(ps);
   }

}


generated quantities {
   real ps_SM[N];
   // This computes respective densities of each mixture component
   for(n in 1:N){
      int k = site[n];
      vector[3] ps;

      ps[1] = prev_comp[k][1]*exp(multi_normal_lpdf(y[n] | mu[1], Sigma[1]));
      ps[2] = prev_comp[k][2]*exp(multi_normal_lpdf(y[n] | mu[2], Sigma[2]));
      ps[3] = prev_comp[k][3]*exp(multi_normal_lpdf(y[n] | mu[3], Sigma[3]));

      // normalise to get P_SM
      ps_SM[n] = ps[3]/sum(ps);
   }
}



