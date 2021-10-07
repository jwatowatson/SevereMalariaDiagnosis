roc_analysis = function(thetas, stan_dat, xvals,
                        xnames, cutofftype, mycuts=NULL){
  Ntest = stan_dat$Ntest
  c_ind = stan_dat$cluster

  sens = sp = xs = mu_vals = sd_vals = list()
  for(i in 1:Ntest){
    xs[[i]] = seq(xvals[i,1], xvals[i,2], length.out = 100)
    mu_vals[[i]] = colMeans(thetas$mu[,i,])
    sd_vals[[i]] = colMeans(thetas$sigma[,i,])

    if(cutofftype[i] == 'upper'){
      sens[[i]] = pnorm(xs[[i]],
                        mu_vals[[i]][c_ind[i,1]],
                        sd_vals[[i]][c_ind[i,1]])
      sp[[i]] = 1-pnorm(xs[[i]],
                        mu_vals[[i]][c_ind[i,2]],
                        sd_vals[[i]][c_ind[i,2]])
      if(!is.null(mycuts)){
        print(xnames[i])
        print(round(100*pnorm(mycuts[i],
                              mu_vals[[i]][c_ind[i,1]],
                              sd_vals[[i]][c_ind[i,1]])))
        print(round(100*(1-pnorm(mycuts[i],
                                 mu_vals[[i]][c_ind[i,2]],
                                 sd_vals[[i]][c_ind[i,2]]))))
      }
    }
    if(cutofftype[i] == 'lower'){
      sens[[i]] = 1-pnorm(xs[[i]],
                          mu_vals[[i]][c_ind[i,1]],
                          sd_vals[[i]][c_ind[i,1]])
      sp[[i]] = pnorm(xs[[i]],
                      mu_vals[[i]][c_ind[i,2]],
                      sd_vals[[i]][c_ind[i,2]])
      if(!is.null(mycuts)){
        print(xnames[i])
        print(round(100*(1-pnorm(mycuts[i],
                                 mu_vals[[i]][c_ind[i,1]],
                                 sd_vals[[i]][c_ind[i,1]]))))
        print(round(100*pnorm(mycuts[i],
                              mu_vals[[i]][c_ind[i,2]],
                              sd_vals[[i]][c_ind[i,2]])))
      }
    }

  }

  for(i in 1:Ntest){
    plot(10^xs[[i]], 100*sens[[i]], type='l',
         panel.first=grid(),log='x',
         xlab=xnames[i],ylab='%',lwd=2, lty=1)
    lines(10^xs[[i]], 100*sp[[i]], type='l',
          log='x', lwd=2, lty=2)
    legend('left', lwd=2, lty=1:2, inset=0.04,
           legend = c('Sensitivity','Specificity'))
  }

  plot(100*(1-sp[[1]]), 100*sens[[1]], type='l',
       xlab = 'False positive rate (%)', ylim=c(0,100),
       ylab = 'True positive rate (%)', xlim=c(0,100),
       lwd=2, lty=1, panel.first=grid())
  for(i in 2:Ntest){
    lines(100*(1-sp[[i]]), 100*sens[[i]], type='l',lty=i,lwd=2)
  }
  lines(0:100,0:100)
  legend('bottomright', legend = xnames,lwd=2,lty=1:3,
         inset=.03, title = 'Biomarker')

}


roc_output = function(thetas, stan_dat, cutofftype, mycuts){

  Ntest = stan_dat$Ntest
  c_ind = stan_dat$cluster
  mu_vals = sd_vals = list()

  specf = sens = array(dim = Ntest)
  for(i in 1:Ntest){

    mu_vals[[i]] = colMeans(thetas$mu[,i,])
    sd_vals[[i]] = colMeans(thetas$sigma[,i,])

    if(cutofftype[i] == 'upper'){
      sens[i] = pnorm(mycuts[i],
                      mu_vals[[i]][c_ind[i,1]],
                      sd_vals[[i]][c_ind[i,1]])
      specf[i] = 1-pnorm(mycuts[i],
                         mu_vals[[i]][c_ind[i,2]],
                         sd_vals[[i]][c_ind[i,2]])

    }
    if(cutofftype[i] == 'lower'){
      sens[i] = 1-pnorm(mycuts[i],
                        mu_vals[[i]][c_ind[i,1]],
                        sd_vals[[i]][c_ind[i,1]])
      specf[i] = pnorm(mycuts[i],
                       mu_vals[[i]][c_ind[i,2]],
                       sd_vals[[i]][c_ind[i,2]])

    }

  }
  sens = prod(sens)
  specf = 1-prod(1-specf)

  return(c(sens, specf))
}
