make_ROC_SM = function(thetas,
                       stan_dat,
                       xvals,
                       Ntest,
                       cutofftype,
                       SM_comp)
{
  sens = sp = xs = mu_vals = sd_vals = list()
  for(i in 1:Ntest){
    xs[[i]] = seq(xvals[i,1], xvals[i,2], length.out = 100)
    mu_vals[[i]] = colMeans(thetas$mu[,,i])
    sd_vals[[i]] = colMeans(thetas$sigma[,,i])
    
    if(cutofftype[i] == 'upper'){
      sens[[i]] = pnorm(xs[[i]],
                        mu_vals[[i]][SM_comp],
                        sd_vals[[i]][SM_comp])
      sp[[i]] = 1-pnorm(xs[[i]],
                        mu_vals[[i]][SM_comp-1],
                        sd_vals[[i]][SM_comp-1])
    }
    if(cutofftype[i] == 'lower'){
      sens[[i]] = 1-pnorm(xs[[i]],
                          mu_vals[[i]][SM_comp],
                          sd_vals[[i]][SM_comp])
      sp[[i]] = pnorm(xs[[i]],
                      mu_vals[[i]][SM_comp-1],
                      sd_vals[[i]][SM_comp-1])
    }
    
  }
  out = list(xs, sens, sp)
  names(out) = c('xs', 'sens', 'sp')
  return(out)
}


make_ROC_SM_tdist = function(thetas,
                             stan_dat,
                             xvals,
                             Ntest,
                             cutofftype,
                             SM_comp)
{
  sens = sp = xs = mu_vals = sd_vals = list()
  for(i in 1:Ntest){
    xs[[i]] = seq(xvals[i,1], xvals[i,2], length.out = 100)
    mu_vals[[i]] = colMeans(thetas$mu[,,i])
    sd_vals[[i]] = colMeans(thetas$sigma[,,i])
    
    if(cutofftype[i] == 'upper'){
      sens[[i]] = pnorm(xs[[i]],
                        mu_vals[[i]][SM_comp],
                        sd_vals[[i]][SM_comp])
      sp[[i]] = 1-pnorm(xs[[i]],
                        mu_vals[[i]][SM_comp-1],
                        sd_vals[[i]][SM_comp-1])
    }
    if(cutofftype[i] == 'lower'){
      sens[[i]] = 1-pnorm(xs[[i]],
                          mu_vals[[i]][SM_comp],
                          sd_vals[[i]][SM_comp])
      sp[[i]] = pnorm(xs[[i]],
                      mu_vals[[i]][SM_comp-1],
                      sd_vals[[i]][SM_comp-1])
    }
    
  }
  out = list(xs, sens, sp)
  names(out) = c('xs', 'sens', 'sp')
  return(out)
}

