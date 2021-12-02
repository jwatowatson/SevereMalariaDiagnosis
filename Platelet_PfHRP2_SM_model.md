---
title: "Platelet count and plasma PfHRP2 concentration to diagnose severe malaria"
author: "James Watson"
date: "5/10/2021"
output: 
  html_document: 
    toc: yes
    number_sections: yes
    keep_md: yes
---





```
## Loading required package: nlme
```

```
## This is mgcv 1.8-36. For overview type 'help("mgcv-package")'.
```

```
## Loading required package: Matrix
```

```
## 
## Attaching package: 'lme4'
```

```
## The following object is masked from 'package:nlme':
## 
##     lmList
```

```
## Package 'mclust' version 5.4.7
## Type 'citation("mclust")' for citing this R package in publications.
```

```
## 
## Attaching package: 'mclust'
```

```
## The following object is masked from 'package:mgcv':
## 
##     mvn
```

```
## Loading required package: StanHeaders
```

```
## Loading required package: ggplot2
```

```
## rstan (Version 2.21.2, GitRev: 2e1f913d3ca3)
```

```
## For execution on a local, multicore CPU with excess RAM we recommend calling
## options(mc.cores = parallel::detectCores()).
## To avoid recompilation of unchanged Stan programs, we recommend calling
## rstan_options(auto_write = TRUE)
```

```
## 
## Attaching package: 'gtools'
```

```
## The following object is masked from 'package:mgcv':
## 
##     scat
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:nlme':
## 
##     collapse
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```
## 
## Attaching package: 'mvtnorm'
```

```
## The following object is masked from 'package:mclust':
## 
##     dmvnorm
```

```
## 
## Attaching package: 'LaplacesDemon'
```

```
## The following objects are masked from 'package:mvtnorm':
## 
##     dmvt, rmvt
```

```
## The following objects are masked from 'package:gtools':
## 
##     ddirichlet, logit, rdirichlet
```

```
## The following objects are masked from 'package:mgcv':
## 
##     dmvn, rmvn
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS  10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] LaplacesDemon_16.1.6 mvtnorm_1.1-2        reshape2_1.4.4      
##  [4] dplyr_1.0.7          gtools_3.9.2         rstan_2.21.2        
##  [7] ggplot2_3.3.5        StanHeaders_2.21.0-7 mclust_5.4.7        
## [10] lme4_1.1-27.1        Matrix_1.3-4         mgcv_1.8-36         
## [13] nlme_3.1-153         RColorBrewer_1.1-2  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.7         lattice_0.20-44    prettyunits_1.1.1  ps_1.6.0          
##  [5] assertthat_0.2.1   digest_0.6.27      utf8_1.2.2         V8_3.4.2          
##  [9] plyr_1.8.6         R6_2.5.1           stats4_4.0.2       evaluate_0.14     
## [13] pillar_1.6.2       rlang_0.4.11       curl_4.3.2         minqa_1.2.4       
## [17] rstudioapi_0.13    callr_3.7.0        nloptr_1.2.2.2     jquerylib_0.1.4   
## [21] rmarkdown_2.11     splines_4.0.2      stringr_1.4.0      loo_2.4.1         
## [25] munsell_0.5.0      compiler_4.0.2     xfun_0.26          pkgconfig_2.0.3   
## [29] pkgbuild_1.2.0     htmltools_0.5.2    tidyselect_1.1.1   tibble_3.1.4      
## [33] gridExtra_2.3      codetools_0.2-18   matrixStats_0.61.0 fansi_0.5.0       
## [37] crayon_1.4.1       withr_2.4.2        MASS_7.3-54        grid_4.0.2        
## [41] jsonlite_1.7.2     gtable_0.3.0       lifecycle_1.0.0    DBI_1.1.1         
## [45] magrittr_2.0.1     scales_1.1.1       RcppParallel_5.1.4 cli_3.0.1         
## [49] stringi_1.7.4      bslib_0.3.0        ellipsis_0.3.2     generics_0.1.0    
## [53] vctrs_0.3.8        boot_1.3-28        tools_4.0.2        glue_1.4.2        
## [57] purrr_0.3.4        parallel_4.0.2     processx_3.5.2     fastmap_1.1.0     
## [61] yaml_2.2.1         inline_0.3.19      colorspace_2.0-2   knitr_1.34        
## [65] sass_0.4.0
```


# Load dataset 

## Check outliers


```
## The dataset contains 2649 patients with measured PfHRP2 and measured platelet counts from 4 studies
```

```
## Patients per study:
```

```
## 
##       Bangladesh   FEAST (Uganda) Kampala (Uganda)   Kilifi (Kenya) 
##              172              567              492             1418
```


## Some data cleaning

![](Platelet_PfHRP2_SM_model_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```
## 
## FALSE  TRUE 
##     6  2643
```

```
## 
##     Bangladesh FEAST (Uganda) Kilifi (Kenya) 
##              1              8             18
```

![](Platelet_PfHRP2_SM_model_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

```
## A total of 27 samples have zero PfHRP2 but more than 1000 parasites per uL
```

```
## After excluding the HRP2 outliers, the dataset contains 2622 patients with measured PfHRP2 and measured platelet counts
```

## Overview of patient characteristics

Results for Table 1 in the paper

```
## 
##       Bangladesh   FEAST (Uganda) Kampala (Uganda)   Kilifi (Kenya) 
##              171              559              492             1400
```

```
##    
##     Bangladesh FEAST (Uganda) Kampala (Uganda) Kilifi (Kenya)
##   0          0            227                0              0
##   1        171            332              492           1400
```

```
##              study age.lower age.median age.upper
## 1       Bangladesh      23.5       30.0      45.0
## 2   FEAST (Uganda)       1.2        2.0       3.3
## 3 Kampala (Uganda)       2.2        3.3       4.6
## 4   Kilifi (Kenya)       1.4        2.4       3.7
```

```
##              study hrp2
## 1       Bangladesh  171
## 2   FEAST (Uganda)  559
## 3 Kampala (Uganda)  492
## 4   Kilifi (Kenya) 1400
```

```
##              study outcome
## 1       Bangladesh    26.9
## 2   FEAST (Uganda)    11.4
## 3 Kampala (Uganda)     6.7
## 4   Kilifi (Kenya)    11.1
```

```
##              study platelet.25% platelet.50% platelet.75%
## 1       Bangladesh         27.0         50.0        139.0
## 2   FEAST (Uganda)         74.5        165.0        326.0
## 3 Kampala (Uganda)         49.0         96.0        169.5
## 4   Kilifi (Kenya)         64.0        111.0        215.0
```

```
##              study  hrp2.25%  hrp2.50%  hrp2.75%
## 1       Bangladesh 1082.9050 2667.0400 6127.5550
## 2   FEAST (Uganda)    0.0000  174.7100 1952.6900
## 3 Kampala (Uganda)  588.0000 1838.4000 4097.4000
## 4   Kilifi (Kenya)  418.7393 2206.7408 5071.5299
```

```
##              study wbc.25% wbc.50% wbc.75%
## 1       Bangladesh   6.900   9.000  11.000
## 2   FEAST (Uganda)   8.400  11.950  18.675
## 3 Kampala (Uganda)   7.500  10.400  15.300
## 4   Kilifi (Kenya)   8.900  12.550  19.000
```

```
##              study para.25% para.50% para.75%
## 1       Bangladesh    23550   148874   348540
## 2   FEAST (Uganda)     3640    37600   153680
## 3 Kampala (Uganda)    10635    42530   198540
## 4   Kilifi (Kenya)     6099    69824   316350
```

```
##                   
##                       1    2    3
##   Bangladesh          0    0    0
##   FEAST (Uganda)    466   46   21
##   Kampala (Uganda)  463    4   23
##   Kilifi (Kenya)   1348   41    7
```

```
##   as.numeric(dat_all$platelet <= 150) hrp2.25% hrp2.50% hrp2.75%
## 1                                   0       24      269     1043
## 2                                   1     1261     3031     6035
```



Correlation between the platelet count and the PfHRP2 concentration


```
## African sites: correlation:
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  log10(dat_all$platelet[ind_Africa]) and log10(dat_all$hrp2 + 1)[ind_Africa]
## t = -31.892, df = 2449, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.5690929 -0.5131183
## sample estimates:
##        cor 
## -0.5417058
```

```
## [1] -186.2477
```

```
##   cor 
## -0.54
```

```
## Bangladesh: correlation:
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  log10(dat_all$platelet[!ind_Africa]) and log10(dat_all$hrp2 + 1)[!ind_Africa]
## t = -4.9404, df = 169, p-value = 1.863e-06
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.4797402 -0.2167256
## sample estimates:
##        cor 
## -0.3552439
```

```
## [1] -5.72988
```

```
##   cor 
## -0.36
```

Summary plot of the biomarker data

![](Platelet_PfHRP2_SM_model_files/figure-html/platelet_hrp2_data-1.png)<!-- -->![](Platelet_PfHRP2_SM_model_files/figure-html/platelet_hrp2_data-2.png)<!-- -->


# Basic data exploration: clustering with mclust

mclust is a generic Bayesian clustering algorithm (fits multivariate normals to the data)

We merge all the data into one and fit mclust


```
## [1] 2477
```

```
## 
##       Bangladesh   FEAST (Uganda) Kampala (Uganda)   Kilifi (Kenya) 
##              170              425              484             1398
```

```
## ---------------------------------------------------- 
## Gaussian finite mixture model fitted by EM algorithm 
## ---------------------------------------------------- 
## 
## Mclust VEV (ellipsoidal, equal shape) model with 9 components: 
## 
##  log-likelihood    n df      BIC       ICL
##       -6408.059 2477 73 -13386.6 -15021.76
## 
## Clustering table:
##   1   2   3   4   5   6   7   8   9 
## 417  27 459 196 230 433 141 305 269
```

```
## ---------------------------------------------------- 
## Gaussian finite mixture model fitted by EM algorithm 
## ---------------------------------------------------- 
## 
## Mclust VVV (ellipsoidal, varying volume, shape, and orientation) model with 4
## components: 
## 
##  log-likelihood    n df       BIC       ICL
##       -3939.675 2622 23 -8060.399 -8886.706
## 
## Clustering table:
##    1    2    3    4 
##  776   62  182 1602
```

![](Platelet_PfHRP2_SM_model_files/figure-html/mclust-1.png)<!-- -->

```
## ---------------------------------------------------- 
## Gaussian finite mixture model fitted by EM algorithm 
## ---------------------------------------------------- 
## 
## Mclust VVV (ellipsoidal, varying volume, shape, and orientation) model with 2
## components: 
## 
##  log-likelihood    n df       BIC       ICL
##       -6667.391 2477 11 -13420.74 -13994.61
## 
## Clustering table:
##    1    2 
##  942 1535
```

![](Platelet_PfHRP2_SM_model_files/figure-html/mclust-2.png)<!-- -->

```
## ---------------------------------------------------- 
## Gaussian finite mixture model fitted by EM algorithm 
## ---------------------------------------------------- 
## 
## Mclust VVV (ellipsoidal, varying volume, shape, and orientation) model with 2
## components: 
## 
##  log-likelihood    n df       BIC       ICL
##       -4785.857 2477 11 -9657.677 -10326.46
## 
## Clustering table:
##    1    2 
## 1148 1329
```

![](Platelet_PfHRP2_SM_model_files/figure-html/mclust-3.png)<!-- -->

```
##                   
##                      1   2   3   4
##   Bangladesh        24   3   0 144
##   FEAST (Uganda)   147  26 179 207
##   Kampala (Uganda) 149   8   0 335
##   Kilifi (Kenya)   456  25   3 916
```


In conclusion, the distribution of the parasite count is less easily decomposed than that of the platelet count or HRP2 concentration. The HRP2 and platelet count is the only one where the estimated break is not orthogonal to either variable.




# Fitting a mixture model to platelets and HRP2
## Two component mixture - not including FEAST

### Main model

Analysis not including the FEAST study - only severe malaria studies

We convert the platelet counts and HRP2 measurements to log10 scale. The Platelet counts then get multiplied by minus 1: this is so that increasing values correspond to more likely severe malaria. This is because the underlying stan model uses the *ordered* vector type to avoid label switching problems. 

It was not possible to get reasonable fits from this model when we included the parasite density. This most likely results from model mis-specification. The parasite density has a lower discriminatory value for severe malaria as shown by our exploratory analyses (see github). Parasite density is mechanistically related to the \textit{Pf}HRP2 concentration. Pre-referral treatment with antimalarials (which is common in severe malaria referrals) can reduce the admission parasite density substantially (thus a much larger variance in the biomarker values for the true severe malaria subgroup). For these reasons, the main models only used the platelet counts and \textit{Pf}HRP2 concentrations.


```
## site_index_SMstudies
##    1    2    3 
##  171  492 1400
```

```
## [1] 2063
```

```
## Priors on mean biomarker values:
```

```
##          [,1]     [,2]
## [1,] -2.30103 -2.00000
## [2,]  2.69897  3.39794
```

```
## Priors on standard deviations around biomarker values:
```

```
##      [,1] [,2]
## [1,]  0.1  0.1
## [2,]  0.1  0.1
```

```
## Priors on prevalence of SM:
```

```
##      [,1] [,2]
## [1,]   19    1
## [2,]   14    6
## [3,]   14    6
```


compile the stan model



Run the stan models - outputs are stored in Rout


We check convergence with the traceplots

![](Platelet_PfHRP2_SM_model_files/figure-html/check_model_fits-1.png)<!-- -->![](Platelet_PfHRP2_SM_model_files/figure-html/check_model_fits-2.png)<!-- -->![](Platelet_PfHRP2_SM_model_files/figure-html/check_model_fits-3.png)<!-- -->![](Platelet_PfHRP2_SM_model_files/figure-html/check_model_fits-4.png)<!-- -->

```
## prevalence estimates for the three studies: model without correlation:
```

```
## [1] 93 69 64
```

```
## prevalence estimates for the three studies: model with correlation:
```

```
## [1] 94 68 59
```

```
## prevalence estimates for the three studies: model with correlation and weak priors:
```

```
## [1] 94 72 64
```

```
## prevalence estimates for the three studies: model with correlation and t-distributions:
```

```
## [1] 95 70 63
```

```
## mean values estimates for the 4 models:
```

```
## thetas_all
##                 
##                  not SM   SM
##   Platelet count    227   70
##   PfHRP2            225 3432
## thetas_all_cor
##                 
##                  not SM   SM
##   Platelet count    197   71
##   PfHRP2            303 3300
## thetas_all_cor_WP
##                 
##                  not SM   SM
##   Platelet count    214   73
##   PfHRP2            211 3252
## thetas_all_cor_tdist
##                 
##                  not SM   SM
##   Platelet count    225   71
##   PfHRP2            281 3275
```

```
## Standard deviation estimates for the 4 models:
```

```
##       
##             [,1]      [,2]
##   [1,] 0.2819436 0.7499934
##   [2,] 0.3032643 0.4416046
```

```
##       
##             [,1]      [,2]
##   [1,] 0.3170179 0.8449156
##   [2,] 0.3170527 0.4391073
```

```
##       
##             [,1]      [,2]
##   [1,] 0.3091195 0.7896354
##   [2,] 0.3182576 0.4512951
```

```
##       
##             [,1]      [,2]
##   [1,] 0.2436313 0.7532043
##   [2,] 0.2863234 0.4214951
```

```
## Correlation in Severe Malaria:
```

```
## [1] 0.3031062
```

```
## Correlation in not Severe Malaria:
```

```
## [1] 0.2345801
```

![](Platelet_PfHRP2_SM_model_files/figure-html/check_model_fits-5.png)<!-- -->

```
## [1] 1077
```

```
## [1] 448
```


Classifications under the bivariate normal model

![](Platelet_PfHRP2_SM_model_files/figure-html/SM_group_classifications_SMonly-1.png)<!-- -->


Classifications under the bivariate t-distribution model

![](Platelet_PfHRP2_SM_model_files/figure-html/SM_group_classifications_SMonly_tdist-1.png)<!-- -->


### Platelet counts in not Severe Malaria and Severe Malaria

![](Platelet_PfHRP2_SM_model_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


### Sensitivity analysis - only African studies

Just use the two severe malaria studies from Africa

```
## site_index_SMAfrica
##    1    2 
##  492 1400
```

Run model and check convergence

```
## 'pars' not specified. Showing first 10 parameters by default.
```

![](Platelet_PfHRP2_SM_model_files/figure-html/fit_stan_model_sensitivityanalysis-1.png)<!-- -->

```
## [1] 65 58
```

```
## Mean estimated values across the groups:
```

```
##                 
##                  not SM   SM
##   Platelet count    195   73
##   PfHRP2            305 3545
```


### ROC curves for main model fit


![](Platelet_PfHRP2_SM_model_files/figure-html/roc_curves_mod2-1.png)<!-- -->

```
## Platelet cutoff of 150,000 per uL has a sensitivity of 85% and a specificity of 65%
```

```
## PLasma PfHRP2 cutoff of 1000 ng/ml has a sensitivity of 88% and a specificity of 73%
```


Empirical ROC curves derived from patient probabilties
![](Platelet_PfHRP2_SM_model_files/figure-html/roc_curves_mod2_empirical-1.png)<!-- -->

### ROC curves for sensitivity analyses

![](Platelet_PfHRP2_SM_model_files/figure-html/sensitivity_roc-1.png)<!-- -->



## Three component mixture - including FEAST


This does not fit when we add parasite counts. So we are only using two biomarkers







```
## site_index_all_studies
##    1    2    3    4 
##  171  559  492 1400
```

```
## Priors on mean biomarker values:
```

```
##           [,1]     [,2]     [,3]
## [1,] -2.477121 -2.30103 -2.00000
## [2,]  0.000000  2.69897  3.39794
```

```
## Priors on standard deviations around biomarker values:
```

```
##      [,1] [,2] [,3]
## [1,]  0.1  0.1  0.1
## [2,]  0.2  0.2  0.2
```

```
## Priors on prevalence of SM:
```

```
##      [,1] [,2] [,3]
## [1,]  0.1    1   19
## [2,]  1.0    6   14
## [3,]  1.0    6   14
## [4,]  3.0    3    3
```

Run models:



![](Platelet_PfHRP2_SM_model_files/figure-html/summary_fit3-1.png)<!-- -->

```
## Estimated prevalences across the 4 studies:
```

```
##        
##         Bangladesh FEAST (Uganda) Kampala (Uganda) Kilifi (Kenya)
##   2.5%        89.9           30.5             65.3           59.2
##   50%         95.6           35.9             72.5           64.7
##   97.5%       98.9           41.1             79.2           69.9
```

```
##                 
##                  not SM not SM   SM
##   Platelet count    286    203   73
##   PfHRP2              1    197 3164
```

```
## In SM the 95% pred interval for the platelet count is 17 to 308 with a mean of 73
```

```
## In SM the 95% pred interval for the PfHRP2 concentration is 414 to 24196 with a mean of 3164
```

```
## In not SM the 95% pred interval for the platelet count is 50 to 832 with a mean of 203
```

```
## In not SM the 95% pred interval for the PfHRP2 concentration is 5 to 7456 with a mean of 197
```




Make a P(SM) variable for subsequent analysis (continuous) and a discrete classification based on a cutoff of 0.5
![](Platelet_PfHRP2_SM_model_files/figure-html/unnamed-chunk-6-1.png)<!-- -->![](Platelet_PfHRP2_SM_model_files/figure-html/unnamed-chunk-6-2.png)<!-- -->



# Downstream analyses

## Severe malaria probability and subphenotypes

Relationship with the two main subphenotypes, cerebral malaria and severe malaria anaemia

```
## False diagnosis rate by coma status and by study:
```

```
##   coma            study        SM
## 1    0       Bangladesh 0.9841270
## 2    1       Bangladesh 0.9722222
## 3    0   FEAST (Uganda) 0.3400901
## 4    1   FEAST (Uganda) 0.4403670
## 5    0 Kampala (Uganda) 0.6447368
## 6    1 Kampala (Uganda) 0.8787879
## 7    0   Kilifi (Kenya) 0.7353952
## 8    1   Kilifi (Kenya) 0.6502525
```

```
## True diagnosis rate by severe anaemia status (<15% haematocrit) and by study:
```

```
##   SMA            study        SM
## 1   0       Bangladesh 0.9753086
## 2   1       Bangladesh 1.0000000
## 3   0   FEAST (Uganda) 0.3454106
## 4   1   FEAST (Uganda) 0.4125874
## 5   0 Kampala (Uganda) 0.8644068
## 6   1 Kampala (Uganda) 0.6835938
## 7   0   Kilifi (Kenya) 0.6469003
## 8   1   Kilifi (Kenya) 0.8596491
```


Main model fits - three component model

```
## , , study = Bangladesh
## 
##       class
## sickle   0   1
##      1   0   0
##      2   0   0
##      3   0   0
## 
## , , study = FEAST (Uganda)
## 
##       class
## sickle   0   1
##      1 281 185
##      2  40   6
##      3  20   1
## 
## , , study = Kampala (Uganda)
## 
##       class
## sickle   0   1
##      1  95 368
##      2   1   3
##      3  17   6
## 
## , , study = Kilifi (Kenya)
## 
##       class
## sickle   0   1
##      1 399 949
##      2  26  15
##      3   6   1
```

![](Platelet_PfHRP2_SM_model_files/figure-html/SM_group_classifications-1.png)<!-- -->

## Parasite densities


```
## [1] "FEAST (Uganda)"
##        1 
## 15.42788
```

```
## [1] "Kampala (Uganda)"
##        1 
## 11.49091
```

```
## [1] "Kilifi (Kenya)"
##       1 
## 13.7813
```

![](Platelet_PfHRP2_SM_model_files/figure-html/parasite_densities-1.png)<!-- -->

```
## [1] "Bangladesh"
##         1 
## 0.4286614
```

```
##                          Estimate Std. Error    t value      Pr(>|t|)
## (Intercept)            3.81224652 0.09267166 41.1371359 1.441362e-277
## P_SM                   1.10526619 0.05676724 19.4701431  2.892819e-78
## studyFEAST (Uganda)   -0.22016532 0.10001509 -2.2013209  2.781257e-02
## studyKampala (Uganda) -0.09041491 0.08812500 -1.0259849  3.050069e-01
## studyKilifi (Kenya)    0.07153498 0.08133515  0.8795088  3.792176e-01
```

```
## The expected fold increase in parasite density between true severe malaria (P_SM=1) and not severe malaria (P_SM=0) is:
```

```
## [1]  9.9 12.7 16.5
```



## Admission haematocrits


![](Platelet_PfHRP2_SM_model_files/figure-html/haematocrits-1.png)<!-- -->

```
##   SM            study   hct
## 1  0       Bangladesh 35.00
## 2  1       Bangladesh 27.35
## 3  0   FEAST (Uganda) 30.00
## 4  1   FEAST (Uganda) 19.65
## 5  0 Kampala (Uganda) 11.60
## 6  1 Kampala (Uganda) 15.50
## 7  0   Kilifi (Kenya) 27.90
## 8  1   Kilifi (Kenya) 19.50
```

```
##   SM            study  hct
## 1  0       Bangladesh  0.0
## 2  1       Bangladesh  4.2
## 3  0   FEAST (Uganda) 22.3
## 4  1   FEAST (Uganda) 26.7
## 5  0 Kampala (Uganda) 71.7
## 6  1 Kampala (Uganda) 46.2
## 7  0   Kilifi (Kenya)  9.2
## 8  1   Kilifi (Kenya) 24.5
```

```
## [1] 19
```

```
## [1] 31
```

```
## [1] 16.4
```

```
## [1] 11.6
```

```
## [1] 19.3
```

![](Platelet_PfHRP2_SM_model_files/figure-html/haematocrits-2.png)<!-- -->

```
## [1] 28.6
```


## Mortality


P(Severe malaria) and mortality

```
## 
## Call:
## glm(formula = outcome ~ P_SM, family = "binomial", data = dat_all[ind, 
##     ])
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -0.5061  -0.5061  -0.4959  -0.4708   2.1277  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  -1.9906     0.1710 -11.639   <2e-16 ***
## P_SM         -0.1640     0.3307  -0.496     0.62    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 397.79  on 558  degrees of freedom
## Residual deviance: 397.54  on 557  degrees of freedom
## AIC: 401.54
## 
## Number of Fisher Scoring iterations: 4
```

```
## 
## Call:
## glm(formula = outcome ~ P_SM, family = "binomial", data = dat_all[ind, 
##     ])
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -0.4682  -0.4505  -0.3993  -0.1943   3.1004  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  -4.8814     0.9188  -5.313 1.08e-07 ***
## P_SM          2.7292     1.0106   2.701  0.00692 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 242.07  on 491  degrees of freedom
## Residual deviance: 229.90  on 490  degrees of freedom
## AIC: 233.9
## 
## Number of Fisher Scoring iterations: 6
```

```
## 
## Call:
## glm(formula = outcome ~ P_SM, family = "binomial", data = dat_all[ind, 
##     ])
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -0.5620  -0.5086  -0.4520  -0.4446   2.1798  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  -1.7657     0.1518 -11.629   <2e-16 ***
## P_SM         -0.5136     0.2134  -2.406   0.0161 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 974.42  on 1399  degrees of freedom
## Residual deviance: 968.77  on 1398  degrees of freedom
## AIC: 972.77
## 
## Number of Fisher Scoring iterations: 4
```

```
## 
## Call:
## glm(formula = outcome ~ P_SM, family = "binomial", data = dat_all[ind, 
##     ])
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -0.8204  -0.8189  -0.8060   1.5830   1.7791  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)
## (Intercept)   -3.224      2.114  -1.525    0.127
## P_SM           2.308      2.168   1.065    0.287
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 199.14  on 170  degrees of freedom
## Residual deviance: 197.40  on 169  degrees of freedom
## AIC: 201.4
## 
## Number of Fisher Scoring iterations: 5
```

![](Platelet_PfHRP2_SM_model_files/figure-html/mortality-1.png)<!-- -->


Summary of lab values by inferred subgroup

```
##   SM            study  age
## 1  0       Bangladesh 38.0
## 2  1       Bangladesh 30.0
## 3  0   FEAST (Uganda)  1.8
## 4  1   FEAST (Uganda)  2.4
## 5  0 Kampala (Uganda)  3.4
## 6  1 Kampala (Uganda)  3.2
## 7  0   Kilifi (Kenya)  2.3
## 8  1   Kilifi (Kenya)  2.4
```

```
##   SM            study log10_parasites
## 1  0       Bangladesh          100000
## 2  1       Bangladesh           79433
## 3  0   FEAST (Uganda)             251
## 4  1   FEAST (Uganda)           39811
## 5  0 Kampala (Uganda)            7943
## 6  1 Kampala (Uganda)           50119
## 7  0   Kilifi (Kenya)           10000
## 8  1   Kilifi (Kenya)           79433
```

```
##   SM            study log10_wbc
## 1  0       Bangladesh         8
## 2  1       Bangladesh         8
## 3  0   FEAST (Uganda)        13
## 4  1   FEAST (Uganda)        10
## 5  0 Kampala (Uganda)        13
## 6  1 Kampala (Uganda)        10
## 7  0   Kilifi (Kenya)        16
## 8  1   Kilifi (Kenya)        13
```

```
##   SM            study sickle.1 sickle.2 sickle.3
## 1  0   FEAST (Uganda)      281       40       20
## 2  1   FEAST (Uganda)      185        6        1
## 3  0 Kampala (Uganda)       95        1       17
## 4  1 Kampala (Uganda)      368        3        6
## 5  0   Kilifi (Kenya)      399       26        6
## 6  1   Kilifi (Kenya)      949       15        1
```

```
##   SM            study outcome
## 1  0       Bangladesh     0.0
## 2  1       Bangladesh    27.5
## 3  0   FEAST (Uganda)    12.0
## 4  1   FEAST (Uganda)    10.4
## 5  0 Kampala (Uganda)     2.7
## 6  1 Kampala (Uganda)     7.9
## 7  0   Kilifi (Kenya)    14.3
## 8  1   Kilifi (Kenya)     9.6
```

## Total white blood cell counts


![](Platelet_PfHRP2_SM_model_files/figure-html/wbcs-1.png)<!-- -->

```
##                           Estimate Std. Error    t value      Pr(>|t|)
## (Intercept)            1.286060464 0.04589606  28.021151 2.645131e-151
## studyFEAST (Uganda)   -0.124664570 0.04376642  -2.848407  4.428221e-03
## studyKampala (Uganda) -0.137776566 0.04229067  -3.257848  1.136974e-03
## studyKilifi (Kenya)   -0.069168701 0.04220474  -1.638885  1.013576e-01
## SM                    -0.114044930 0.01120286 -10.179985  6.756193e-24
## age                   -0.006774506 0.00118929  -5.696263  1.361504e-08
```

```
## [1] 1.300304
```

```
## WBC > 15,000 for SM versus not SM:
```

```
## 
## Family: binomial 
## Link function: logit 
## 
## Formula:
## as.numeric(wbc > 15) ~ study + SM + s(age, k = 3)
## 
## Parametric coefficients:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)            2.02403    0.46847   4.320 1.56e-05 ***
## studyFEAST (Uganda)   -2.75440    0.51248  -5.375 7.68e-08 ***
## studyKampala (Uganda) -2.71246    0.48982  -5.538 3.07e-08 ***
## studyKilifi (Kenya)   -2.36279    0.49738  -4.751 2.03e-06 ***
## SM                    -0.71074    0.09796  -7.256 3.99e-13 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##          edf Ref.df Chi.sq p-value    
## s(age) 1.976  1.999  70.62  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0775   Deviance explained = 6.37%
## UBRE = 0.18848  Scale est. = 1         n = 2393
```

```
## Odds-ratio:
```

```
## [1] 0.45 0.49 0.54
```

```
##                   
##                    FALSE TRUE
##   Bangladesh         149   22
##   FEAST (Uganda)     354  204
##   Kampala (Uganda)   365  126
##   Kilifi (Kenya)     870  530
```

## Bacteraemia

Relationship between the blood culture data and the probability of severe malaria.
We also look at whether WBC values predict positive blood culture


```
## Blood culture positive numbers by study:
```

```
## , ,  = Bangladesh
## 
##    
##       0   1
##   0   0   0
##   1   0   0
## 
## , ,  = FEAST (Uganda)
## 
##    
##       0   1
##   0 273  39
##   1 165  18
## 
## , ,  = Kampala (Uganda)
## 
##    
##       0   1
##   0   0   0
##   1   0   0
## 
## , ,  = Kilifi (Kenya)
## 
##    
##       0   1
##   0 406  29
##   1 943  22
```

```
## Blood culture positive numbers in FEAST (all):
```

```
## 
##    0    1 <NA> 
##  438   57   64
```

```
## Blood culture positive numbers in FEAST (malaria RDT+):
```

```
## 
##    0    1 <NA> 
##  263   35   34
```

```
## Blood culture positive numbers in Kenya:
```

```
## 
##    0    1 
## 1349   51
```

```
## [1] 0.8855098
```

```
## main model for bacteraemia:
```

```
## 
## Call:
## glm(formula = BS ~ SM + study, family = "binomial", data = dat_all[ind_mal, 
##     ])
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -0.6240  -0.3544  -0.2269  -0.2269   2.7102  
## 
## Coefficients:
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)          -1.5375     0.2074  -7.415 1.22e-13 ***
## SM                   -0.9108     0.2257  -4.036 5.44e-05 ***
## studyKilifi (Kenya)  -1.1986     0.2318  -5.171 2.32e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 680.62  on 1697  degrees of freedom
## Residual deviance: 637.41  on 1695  degrees of freedom
##   (697 observations deleted due to missingness)
## AIC: 643.41
## 
## Number of Fisher Scoring iterations: 6
```

```
##               Estimate Std. Error   z value     Pr(>|z|)
## (Intercept) -1.6945957  0.2564364 -6.608249 3.888921e-11
## SM          -0.5781364  0.3614468 -1.599506 1.097081e-01
```

```
##              Estimate Std. Error    z value     Pr(>|z|)
## (Intercept) -2.639057  0.1922129 -13.729863 6.726034e-43
## SM          -1.118966  0.2888922  -3.873301 1.073710e-04
```

```
## Odds ratio for blood culture+ in SM versus not SM:
```

```
## [1] 0.26 0.40 0.63
```

```
## 
## Call:
## glm(formula = BS ~ log10_wbc + study, family = "binomial", data = dat_all[ind_mal, 
##     ])
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -0.5252  -0.2777  -0.2727  -0.2696   2.6240  
## 
## Coefficients:
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)          -2.1515     0.4822  -4.462 8.11e-06 ***
## log10_wbc             0.1247     0.4127   0.302    0.762    
## studyKilifi (Kenya)  -1.2638     0.2303  -5.488 4.08e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 680.62  on 1697  degrees of freedom
## Residual deviance: 653.53  on 1695  degrees of freedom
##   (697 observations deleted due to missingness)
## AIC: 659.53
## 
## Number of Fisher Scoring iterations: 6
```

# Validation using sickle trait


```
## 
## Call:
## glm(formula = SM ~ sickle + study, family = "binomial", data = data.frame(SM = dat_all$SM[ind], 
##     sickle = as.numeric(dat_all$sickle == 2)[ind], study = as.factor(dat_all$study)[ind]))
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.7857  -1.0038   0.6737   0.8384   1.9731  
## 
## Coefficients:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)           -0.42308    0.09322  -4.539 5.66e-06 ***
## sickle                -1.36948    0.25277  -5.418 6.03e-08 ***
## studyKampala (Uganda)  1.79061    0.14762  12.130  < 2e-16 ***
## studyKilifi (Kenya)    1.28798    0.10953  11.759  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 3082.3  on 2367  degrees of freedom
## Residual deviance: 2828.8  on 2364  degrees of freedom
## AIC: 2836.8
## 
## Number of Fisher Scoring iterations: 4
```

```
## Odds ratio for SM for HbAS versus HbAA:
```

```
## [1] 0.15 0.25 0.42
```

```
## [1] "Kilifi (Kenya)"
## [1] 0.13 0.24 0.46
## [1] "Kampala (Uganda)"
## [1] 0.08 0.77 7.53
## [1] "FEAST (Uganda)"
## [1] 0.09 0.23 0.55
```

```
## , ,  = Bangladesh
## 
##    
##       1   2   3
##   0   0   0   0
##   1   0   0   0
## 
## , ,  = FEAST (Uganda)
## 
##    
##       1   2   3
##   0 281  40  20
##   1 185   6   1
## 
## , ,  = Kampala (Uganda)
## 
##    
##       1   2   3
##   0  95   1  17
##   1 368   3   6
## 
## , ,  = Kilifi (Kenya)
## 
##    
##       1   2   3
##   0 399  26   6
##   1 949  15   1
```

```
## 
## Call:
## glm(formula = SM ~ sickle + study, family = "binomial", data = data.frame(SM = dat_all$SM[ind], 
##     sickle = as.numeric(dat_all$sickle == 3)[ind], study = as.factor(dat_all$study)[ind]))
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.7815  -1.0054   0.6765   0.8380   2.4244  
## 
## Coefficients:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)           -0.41900    0.09432  -4.442 8.90e-06 ***
## sickle                -2.46542    0.40554  -6.079 1.21e-09 ***
## studyKampala (Uganda)  1.77702    0.14754  12.045  < 2e-16 ***
## studyKilifi (Kenya)    1.28482    0.11151  11.522  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 3018.5  on 2327  degrees of freedom
## Residual deviance: 2773.9  on 2324  degrees of freedom
## AIC: 2781.9
## 
## Number of Fisher Scoring iterations: 4
```

```
## Odds ratio for SM for HbSS versus HbAA:
```

```
## [1] 0.04 0.08 0.19
```

```
## [1] "Kilifi (Kenya)"
## [1] 0.01 0.07 0.58
## [1] "Kampala (Uganda)"
## [1] 0.03 0.09 0.24
## [1] "FEAST (Uganda)"
## [1] 0.01 0.08 0.57
```

```
## 
## Call:
## glm(formula = SM ~ sickle + study, family = "binomial", data = data.frame(SM = dat_all$SM, 
##     sickle = dat_all$sickle, study = as.factor(dat_all$study)))
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.7894  -1.0023   0.6713   0.8394   2.4717  
## 
## Coefficients:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)             0.8628     0.1951   4.423 9.72e-06 ***
## sickle                 -1.2897     0.1613  -7.994 1.31e-15 ***
## studyKampala (Uganda)   1.8025     0.1452  12.416  < 2e-16 ***
## studyKilifi (Kenya)     1.2889     0.1091  11.816  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 3177.2  on 2418  degrees of freedom
## Residual deviance: 2869.2  on 2415  degrees of freedom
##   (203 observations deleted due to missingness)
## AIC: 2877.2
## 
## Number of Fisher Scoring iterations: 4
```

```
## Odds ratio for SM for additional S allele (additive model):
```

```
## [1] 0.20 0.28 0.38
```




# Compare with previously published probabilities

![](Platelet_PfHRP2_SM_model_files/figure-html/comparison_old_model-1.png)<!-- -->![](Platelet_PfHRP2_SM_model_files/figure-html/comparison_old_model-2.png)<!-- -->

```
## 
## 	Pearson's product-moment correlation
## 
## data:  model1_probs$P_SM1 and dat_Kilifi$P_SM
## t = 31.253, df = 1398, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.6094134 0.6711674
## sample estimates:
##       cor 
## 0.6413278
```





Save output




# Optimal combination of thresholds





Plot results

```
## [1]  5.18 40.44
```

```
## [1] 50.11 92.21
```

```
## [1] 77 94
```

```
## There are 0 combinations that have a true positive rate greater than 75 and a false positive rate less than 10
```

```
##       plt_thresh hrp2_thresh  x1_trans x2_trans    TP    FP ProbSM
## 18231        150        1000 -2.176091        3 75.96 13.43     89
```

![](Platelet_PfHRP2_SM_model_files/figure-html/optimal_combination-1.png)<!-- -->


# Predicting Severe Malaria using HRP2, parasite counts and HCT

Slightly random: predict SM status using HRP2 and HCT

```
## 
## Family: binomial 
## Link function: logit 
## 
## Formula:
## SM ~ s(hct, k = 3) + s(log10(hrp2 + 1), k = 3) + s(log10_parasites, 
##     k = 3) + s(study, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error z value Pr(>|z|)
## (Intercept)   -1.010      6.051  -0.167    0.867
## 
## Approximate significance of smooth terms:
##                      edf Ref.df Chi.sq  p-value    
## s(hct)             1.949  1.997  16.70 0.000224 ***
## s(log10(hrp2 + 1)) 1.978  1.999 290.23  < 2e-16 ***
## s(log10_parasites) 1.747  1.936  16.64 0.000142 ***
## s(study)           2.979  3.000  47.87  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.864   Deviance explained = 83.1%
## UBRE = -0.77747  Scale est. = 1         n = 2474
```

