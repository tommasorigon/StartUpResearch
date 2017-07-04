# Model 1
Alessia Caponera, Francesco Denti, Tommaso Rigon, Andrea Sottosanti and Alan Gelfand  



In this document, we estimated the `model1`, according to our terminology. First of all, we load the necessary data in memory. The data matrix `B` is composed by

- The observations belonging to the first hemisphere (regions from 1 to 35).
- The `10`-th subject.
- The first scan, that is, `k=1`.


```r
rm(list=ls())

# This load the dataset "Y"
load("../Dataset/fMRI-ROI-time-series.RData") 
# Location and names of the Desikan atlas
ROI  <- read.table("../Dataset/ROI-covariates.txt",header=TRUE)
# Information on subjects
SUBJ <- read.table("../Dataset/SUBJ-covariates.txt",header=TRUE)

# This is the selected dataset
B <- Y[1:35,,10,1]
```

The matrix `D` computed below contains the euclidean distances among centroids. Again, we retained only the first `35` regions. In the `W` matrix we stored instead the **proximity matrix** for the CAR prior. The parameters `k1` and `k2` are selected so that the maximum off-diagonal values are between `0.1` and `1`.


```r
D  <- as.matrix(dist(subset(ROI, select=c(centroid_x,centroid_y,centroid_z))))[1:nrow(B),1:nrow(B)]

index_D   <-  matrix(TRUE,nrow(D),nrow(D)); diag(index_D)<- FALSE
d_lower   <- min(D[index_D])
d_upper   <- max(D[index_D])

k1        <- 10^(d_lower/(d_upper - d_lower))
k2        <- log(k1)/d_lower
W         <- k1*exp(-k2*D); diag(W) <- 0
```

The prior hyperparameters are specified through the following list


```r
prior <- list(psi=3/100,              # Scale parameter for the GP
              sigma2_0 = 10^4,        # Variance of the intercept
              a_sigma2=1, b_sigma2=1, # Inverse-Gamma hyperparameters for sigma2
              a_gamma2=1, b_gamma2=1, # Inverse-Gamma hyperparameters for gamma2
              a_tau2=1, b_tau2=1      # Inverse-Gamma hyperparameters for tau2
              )
```


## Gibbs sampling

In the following, the MCMC chain is computed via Gibbs sampling. The required function [`model1.R`](https://github.com/tommasorigon/StartUpResearch/blob/master/Models/model1.R) are made avaible. We remark that it could be quite time consuming. We computed a total of `25000` replications after a period of `500` as burn in. The chain was thinned every `5` replicates for a total of `R=5000` replicates that were retained for inference.

We estimated the model both in-sample and using a fraction of data (about the 75%), in which we randomly eliminated observations.


```r
source("model1.R")

set.seed(123)
# Attention! This part will require a lot of time
fit_model1  <- Gibbs_model1(B=B, W=W, prior=prior, 
                            R=5000, burn_in=500, thinning=5, verbose=TRUE)

# Remove the 25% of values 
index <- sample(nrow(B)*ncol(B),floor(0.25*length(B)))
B_NA <- B; B_NA[index] <- NA 

fit_model1_NA  <- Gibbs_model1(B=B_NA, W=W, prior=prior, 
                            R=5000, burn_in=500, thinning=5, verbose=TRUE)
save.image("model1.RData")
```


## Diagnostic graphs: autocorrelations and traceplots

We displayed some traceplots for the variance parameters in order to assess the mixing and the convergence of the chain.

```r
par(mfrow=c(2,3))

# Autocorrelations
plot(fit_model1$gamma2,type="l")
plot(fit_model1$tau2,type="l")
plot(fit_model1$sigma2,type="l")

# Autocorrelations of hyperparameters
acf(fit_model1$gamma2)
acf(fit_model1$tau2)
acf(fit_model1$sigma2)
```

![](https://github.com/tommasorigon/StartUpResearch/blob/master/Models/Model1_files/figure-html/acf%20and%20traceplots-1.png?raw=true)<!-- -->

Due to the large amount of parameters, the convergence for the spatial random effects and the Gaussian process is monitored via effective sample sizes. Notice that the we computed the ESS of the Gaussian Process **together** with the intercept term.


```r
library(coda)

# Summary of effective sample size for phi
quantile(effectiveSize(as.mcmc(fit_model1$phi)),c(0.25,0.5,0.75))
```

```
##  25%  50%  75% 
## 5000 5000 5000
```

```r
# Summary of effective sample size for the GP + intercept
quantile(effectiveSize(as.mcmc(fit_model1$beta_0 + fit_model1$Z)),c(0.25,0.5,0.75))
```

```
##  25%  50%  75% 
## 5000 5000 5000
```

## Posterior analysis

Below, we report the posterior mean for the spatial coefficients with a credible interval of level 95%.


```r
alpha <- 0.05
object <- fit_model1$phi
tab <-cbind(cbind(apply(object,2,mean),apply(object,2,function(x) quantile(x,alpha/2)),apply(object,2,function(x) quantile(x,1-alpha/2))))

colnames(tab) <- c("Mean","Lower","Upper")
rownames(tab) <- rownames(Y)[1:nrow(B)] # Number of the breain region

knitr::kable(round(t(tab),digits=2),format="markdown")
```



|      | lh-unknown| lh-bankssts| lh-caudalanteriorcingulate| lh-caudalmiddlefrontal| lh-corpuscallosum| lh-cuneus| lh-entorhinal| lh-fusiform| lh-inferiorparietal| lh-inferiortemporal| lh-isthmuscingulate| lh-lateraloccipital| lh-lateralorbitofrontal| lh-lingual| lh-medialorbitofrontal| lh-middletemporal| lh-parahippocampal| lh-paracentral| lh-parsopercularis| lh-parsorbitalis| lh-parstriangularis| lh-pericalcarine| lh-postcentral| lh-posteriorcingulate| lh-precentral| lh-precuneus| lh-rostralanteriorcingulate| lh-rostralmiddlefrontal| lh-superiorfrontal| lh-superiorparietal| lh-superiortemporal| lh-supramarginal| lh-frontalpole| lh-temporalpole| lh-transversetemporal|
|:-----|----------:|-----------:|--------------------------:|----------------------:|-----------------:|---------:|-------------:|-----------:|-------------------:|-------------------:|-------------------:|-------------------:|-----------------------:|----------:|----------------------:|-----------------:|------------------:|--------------:|------------------:|----------------:|-------------------:|----------------:|--------------:|---------------------:|-------------:|------------:|---------------------------:|-----------------------:|------------------:|-------------------:|-------------------:|----------------:|--------------:|---------------:|---------------------:|
|Mean  |        0.0|         0.0|                       0.00|                    0.0|               0.0|      0.00|          0.00|        0.00|                0.00|                0.00|                0.00|                0.00|                    0.00|       0.00|                   0.00|              0.00|               0.00|           0.00|                0.0|             0.00|                0.00|             0.00|            0.0|                   0.0|           0.0|         0.00|                        0.00|                    0.00|               0.00|                0.00|                 0.0|              0.0|           0.00|            0.00|                   0.0|
|Lower |       -0.1|        -0.1|                      -0.11|                   -0.1|              -0.1|     -0.11|         -0.11|       -0.11|               -0.11|               -0.11|               -0.10|               -0.11|                   -0.11|      -0.11|                  -0.11|             -0.10|              -0.10|          -0.11|               -0.1|            -0.11|               -0.11|            -0.11|           -0.1|                  -0.1|          -0.1|        -0.11|                       -0.11|                   -0.11|              -0.11|               -0.11|                -0.1|             -0.1|          -0.12|           -0.11|                  -0.1|
|Upper |        0.1|         0.1|                       0.11|                    0.1|               0.1|      0.12|          0.11|        0.11|                0.11|                0.10|                0.11|                0.12|                    0.11|       0.11|                   0.11|              0.11|               0.11|           0.11|                0.1|             0.11|                0.11|             0.12|            0.1|                   0.1|           0.1|         0.11|                        0.11|                    0.11|               0.11|                0.11|                 0.1|              0.1|           0.12|            0.11|                   0.1|

This table substantially reveals that **no additive spatial structure** seems to be present, thus discouraging the exploration of such a specification. For this reason, we will made available further specification.

The path of the Gaussian process + the intercept is displayed below and compared with the average path.


```r
par(mfrow=c(1,1))
plot(1:404,colMeans(B),type="l")
lines(1:404,colMeans(fit_model1$beta + fit_model1$Z),col="red",lty="dashed")
```

![](https://github.com/tommasorigon/StartUpResearch/blob/master/Models/Model1_files/figure-html/GP%20path-1.png?raw=true)<!-- -->

## Predictive performance

In the following, in-sample and out-of-sample predictive performance are evaluated below.

```r
# IN-SAMPLE predictive performance
pred_model1 <- matrix(0,nrow(B),ncol(B))

for(l in 1:nrow(B)){
  pred_model1[l,] <- colMeans(fit_model1$beta_0 + fit_model1$phi[,l] + fit_model1$Z)
}

RMSE_model1 <- sqrt(mean( (pred_model1 - B)^2))

# OUT-OF-SAMPLE predictive performance
for(l in 1:nrow(B)){
  pred_model1[l,] <- colMeans(fit_model1_NA$beta_0 + fit_model1_NA$phi[,l] + fit_model1_NA$Z)
}

RMSE_model1_NA <- sqrt(mean( (pred_model1[index] - B[index])^2))

tab <- cbind(RMSE_model1,RMSE_model1_NA,sd(c(B)))
colnames(tab) <- c("In-sample RMSE","Out-of-sample RMSE","Standard deviation of the data")
knitr::kable(tab,format="markdown")
```



| In-sample RMSE| Out-of-sample RMSE| Standard deviation of the data|
|--------------:|------------------:|------------------------------:|
|       2.631651|           2.724732|                       3.432074|

As expected, the in-sample RMSE is lower than the out-of-sample RMSE. However, they are both lower than the standard deviation.


