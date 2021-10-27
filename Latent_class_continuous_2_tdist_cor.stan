// This model assumes two normal distributions (one for disease/one for not disease)
// Component 1: not SM
// Component 2: SM
data {
   int N;                                // number of patients
   int Ntest;                            // number of tests
   vector[Ntest] y[N];                   // data
   vector[2] prior_mu[Ntest];
   vector<lower=0>[2] prior_mu_sigma[Ntest];
   int<lower=0> K_sites;                 // number of different sites (populations)
   int<lower=1,upper=K_sites> site[N];   // index of sites for each patient
   vector[2] beta_prior_prev[K_sites];
   real<lower=0> nu;
   int<lower=0> df;
}
parameters {
   real<lower=0,upper=1> prev[K_sites];   //mixing proportions
   ordered[2] mu_test[Ntest];
   vector<lower=0>[Ntest] sigma[2];
   cholesky_factor_corr[Ntest] Lcorr[2];// cholesky factor (L_u matrix for R)
}
transformed parameters{
   corr_matrix[Ntest] R[2]; // correlation matrix
   cov_matrix[Ntest] Sigma[2]; // VCV matrix
   vector[Ntest] mu[2];

   for(t in 1:2){
      for(i in 1:Ntest){
         mu[t][i] = mu_test[i][t];
      }
   }

   for(i in 1:2){
      R[i] = multiply_lower_tri_self_transpose(Lcorr[i]);
      Sigma[i] = quad_form_diag(R[i], sigma[i]);
   }
}

model {
   // **** Prior **** //
   for(t in 1:2){
      for(i in 1:Ntest){
         sigma[t][i] ~ normal(.5, .5);
      }
      Lcorr[t] ~ lkj_corr_cholesky(nu); // prior for cholesky factor
   }
   for(ks in 1:K_sites){
      prev[ks] ~ beta(beta_prior_prev[ks][1],beta_prior_prev[ks][2]);
   }
   for(i in 1:Ntest){
      mu_test[i] ~ normal(prior_mu[i], prior_mu_sigma[i]);
   }
   // **** Likelihood **** //
   for(n in 1:N){
      int k = site[n];
      // group 1 is notSM; group 2 is SM
      target += log_sum_exp(log1m(prev[k])+
      multi_student_t_lpdf(y[n] | df, mu[1], Sigma[1]),
      log(prev[k])+
      multi_student_t_lpdf(y[n] | df, mu[2], Sigma[2]));
   }

}


generated quantities {
   real ps_SM[N];
   // This computes respective densities of each mixture component
   for(n in 1:N){
      int k = site[n];
      real p1;
      real p2;
      p1 = (1-prev[k])*exp(multi_student_t_lpdf(y[n] | df, mu[1], Sigma[1]));
      p2 = prev[k]*exp(multi_student_t_lpdf(y[n] | df, mu[2], Sigma[2]));

      // normalise to get P_SM
      ps_SM[n] = p2/(p1+p2);
   }
}



