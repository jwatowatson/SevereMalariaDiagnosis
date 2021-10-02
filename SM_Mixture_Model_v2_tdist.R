# make a hierarchical structure for theta

require(rstan)
SM_mixture_model_v2_tdist = '
   data {
      int D;                           // dimension of y
      int K;                           // number of gaussians
      int N;                           // number of data
      vector[D] y[N];                  // data
      int<lower=0> K_sites;
      int<lower=1,upper=K_sites> site[N];
      vector[K] prior_mu_platelets;    // prior on mixture mean values
      vector[K] prior_mu_pfhrp2;       // prior on mixture mean values
      vector[K] prior_sd_mu_platelets; // prior on mixture mean values
      vector[K] prior_sd_mu_pfhrp2;    // prior on mixture mean values
      vector[K] alpha_prior[K_sites];
   }
   parameters {
      simplex[K] theta[K_sites];          //mixing proportions
      ordered[K] mu_platelets;            //mixture component means
      ordered[K] mu_pfhrp2;               //mixture component means
      vector<lower=0>[K] sigma_platelets;
      vector<lower=0>[K] sigma_pfhrp2;
      ordered[K] mu_platelets_site[K_sites];
      ordered[K] mu_pfhrp2_site[K_sites];
      cholesky_factor_corr[2] L[K];       //cholesky factor of covariance
      vector<lower=0>[2] Sigma_scale[K];
   }
   transformed parameters {
      matrix[D, D] Sigma[K];
      matrix[D,K] mu[K_sites];
      // covariance matrices
      for(ks in 1:K_sites){
         for(k in 1:K){
            mu[ks][1,k] = mu_platelets_site[ks][k];
            mu[ks][2,k] = mu_pfhrp2_site[ks][K-k+1];
         }
      }
      for(k in 1:K){
         Sigma[k] = crossprod(diag_pre_multiply(Sigma_scale[k], L[k]));
      }
   }
   model {
      real ps[K];
      mu_platelets ~ normal(prior_mu_platelets,prior_sd_mu_platelets);
      mu_pfhrp2 ~ normal(prior_mu_pfhrp2,prior_sd_mu_pfhrp2);
      sigma_platelets ~ normal(.01,.01);
      sigma_pfhrp2 ~ normal(.01,.01);

      for(ks in 1:K_sites){
         theta[ks] ~ dirichlet(alpha_prior[ks]);
         mu_platelets_site[ks] ~ normal(mu_platelets, sigma_platelets);
         mu_pfhrp2_site[ks] ~ normal(mu_pfhrp2, sigma_pfhrp2);
      }

      for(k in 1:K){
         L[k] ~ lkj_corr_cholesky(1);
         Sigma_scale[k] ~ normal(.5, 1);
      }

      for(n in 1:N){
         int ks = site[n];
         for (k in 1:K){
            ps[k] = log(theta[ks][k]) + multi_normal_lpdf(y[n] | mu[ks][,k], Sigma[k]);
         }
         target += log_sum_exp(ps);
      }

   }
   generated quantities {
      vector[K] ps_post[N];
      // This computes respective densities of each mixture component
      for(n in 1:N){
         int ks = site[n];
         vector[K] probs_raw;
         for(k in 1:K){
            probs_raw[k] = theta[ks][k]*exp(multi_normal_lpdf(y[n] | mu[ks][,k], Sigma[k]));
         }
         // normalise to 1
         ps_post[n] = probs_raw/sum(probs_raw);
      }
   }

'


my_mod = rstan::stan_model(model_code = SM_mixture_model_v2_tdist)


