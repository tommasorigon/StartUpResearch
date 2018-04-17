---
title: " Hierarchical spatio-temporal modeling of resting state fMRI: an Application to Multivatiate Time Series Data"
author: "Caponera, Denti, Gelfand, Rigon, Sottosanti"
date: "Last update: April 2018"
output: 
 html_document:
    keep_md: true
---



## Introduction
In this repository we apply the statistical methodology developed in the paper ``Hierarchical spatio-temporal modeling of resting state fMRI'' to a financial dataset, provided in the R package \texttt{MTS}. 
Motivated by \emph{rs-FMRI} data, our work provides a flexible and parsimonious setting for areally-referenced time series. Since in various applications a neighborhood matrix is not available and the idea of which distance to use is not always completely clear, our model does not rely on any distance concept. This assumption makes our model robust with respect to the definition of distance, and hence allows us to apply our methodology to different types of data. One simple example of the  structures that our model is able to handle is given by \emph{multivariate time series}.
More precisely, we use our model to describe the relationships among the time series, providing as output an estimate of the covariance matrix which could be useful for the sake of interpretability and description of the phenomenon under study.

## Data Analysis
To testify the adaptability of model, we apply our setting to describe the covaraince structure among ten time series, referring to simple monthly returns of ten U.S. stocks. As already stated,  The dataset `tenstocks` is available in the package `MTS`.
First of all, since our model was developed for zero-mean time series, we center the data. Let us have a glimpse of our data with the following plot.


```r
# This load the dataset
data("mts-examples")
dd                  <- as.character(tenstocks[,1])
dd                  <- as.Date(dd,"%Y %m %d")
DAT                 <- as.data.frame(apply(tenstocks[,-1],2,function(x) scale(x, center = T, scale = F)))
DAT$date            <- dd
data.plot           <- melt(DAT,id.vars = "date")
colnames(data.plot) <- c("Time","Stock","Value")
data.plot           <- as.tibble(data.plot)
ggplot(data=data.plot, aes(x=Time,y=Value,group=Stock)) + geom_line(alpha=0.60,aes(col=Stock)) + theme_bw() + xlab("Time") + ylab("Value") +ggtitle("Centered monthly simple returns of ten U.S. stocks")
```

![](tenStock_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

One can naively describe the correlation structure computing the Pearson Correlation Coefficient. An example is reported in the following raster plot. 
Spurious correlation, arised by time dependency, could be present in this first result. However, this matrix could be a useful benchmark. 
To support our last claim, we report also the ACF plots (in a time series form, to improve interpretability). 


```r
DAT <- DAT[,-11]

# Pearson correlation
Pear_10 <- cor((DAT))
ggcorrplot(Pear_10,ggtheme=ggplot2::theme_dark,colors=c("red","white","blue"),legend.title="Correlation") +ggtitle("Pearson correlation index among monthly simple returns of ten U.S. stocks")
```

![](tenStock_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
# Autocorrelation
ACF_10           <- matrix(0,21,10)
for(l in 1:10) {ACF_10[,l] <- acf(DAT[,l],plot=FALSE,lag.max=20)$acf[,,1]}
ACF_10           <- as.data.frame(ACF_10)
colnames(ACF_10) <- colnames((DAT))
ACF              <- rbind(melt(ACF_10))
ACF$ind          <- rep(0:20,10)

ggplot(data=ACF,aes(x=ind,y=value,group=variable)) + geom_line(aes(col=variable), alpha=.8,size=0.5) + geom_point(aes(col=variable),size=0.4)+ geom_hline(yintercept=0) + ylab("Autocorrelation") + xlab("Lag") + theme_bw()  + scale_x_continuous(breaks = round(seq(0, 19, by = 1))) + geom_hline(yintercept =qnorm(c(0.025, 0.975))/sqrt(nrow(DAT)),col=2)+ggtitle("Autocorrelation Plot, Maximum lag set to 20")
```

![](tenStock_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

From the correlation matrix, two ``clusters'' of highly correlated time series appear. Moreover, even if it is not so remarkable, some kind of autocorrelation seems to be present among the ten stocks. We then prepare the data for being processed by our model and we import our functions. 


```r
source("/home/fra/Downloads/StartUpResearch-master/TenStockApplication/functions.R")
DAT   <- t(DAT)
n_t   <- ncol(DAT)             # Number of time observations
n_l   <- nrow(DAT)             # Number of generated processes
colnames(DAT) <- 1:n_t;        # Set the names of the colums
time  <- seq(0,ncol(DAT),length=ncol(DAT))
```

In order to perform in-sample prediction we create two mutually exclusive grids of points. The `DAT' dataset is then splitted into two parts: training and test.


```r
set.seed(123)
time_grid <- sort(sample(n_t,round(0.75*n_t))) # The 75% of the time columns are used.
new_grid  <- setdiff(1:n_t,time_grid)          # Grid for prediction
DAT_train <- DAT[,time_grid]
DAT_test  <- DAT[,new_grid]
```

After having splitted the data into train and test dataset, we firstly run 15000 iterations and we discard 5000 of them as burn-in period (thinning step of 5 observations) for $K\in\{2,3,4\}$, as suggested by the previous exploratory plots.


```r
r <- foreach(k=2:3) %dopar% {
  prior             <- list(a_sigma = 0.001, b_sigma = 0.001, gamma2 = 100, K = k, psi= 0.03)
  Sigma_Metropolis  <- cov(bootPPCA(t(DAT_train),prior$K,R=1000))
  tuning            <- list(cholSigmaMetropolis = chol(Sigma_Metropolis))
  Metropolis_MCMC(DAT_train, time_grid, prior, 
                                       Iter=10000, burn_in=5000, start= NULL, 
                                       thinning=5, tuning, verbose = TRUE)
}
#save(r,file="Parallel_results2_3.Rdata")

prior             <- list(a_sigma = 0.001, b_sigma = 0.001, gamma2 = 100, K = 2, psi= 0.03)
Sigma_Metropolis  <- cov(bootPPCA(t(DAT_train),prior$K,R=1000))
tuning            <- list(cholSigmaMetropolis = chol(Sigma_Metropolis))
r2 <- Metropolis_MCMC(DAT_train, time_grid, prior, 
                                       Iter=1000000, burn_in=250000, start= NULL, 
                                       thinning=5, tuning, verbose = TRUE)
#save(r2,file="Parallel_results2.Rdata")
```




```r
# -----------------------------------
# Model summary
# -----------------------------------
model_summary                  <- data.frame(K=c(2:4),rbind(t(unlist(sapply(r,function(x) IC(x)))))[1:3,])
colnames(model_summary)[2:4]   <- c("DIC",  "p","p_DIC")
model_summary$Acceptance_ratio <- c(sapply(r, function(x) x$acceptance_ratio)[1:3])
model_summary$MAP              <- c(sapply(r, function(x) max(x$logpost))[1:3])
model_summary$RMSE_train       <- c(
    sapply(r,function(x) sqrt(mean((DAT_train - predict(x, DAT_train, 1:n_t)[,time_grid])^2))))[1:3]

model_summary$RMSE_test        <- c(  sapply(r,function(x) sqrt(mean((DAT_test - predict(x, DAT_train, 1:n_t)[,new_grid])^2))))[1:3]
knitr::kable(model_summary,digits=4,format = "markdown")
```



|  K|       DIC|  p|   p_DIC| Acceptance_ratio|      MAP| RMSE_train| RMSE_test|
|--:|---------:|--:|-------:|----------------:|--------:|----------:|---------:|
|  2| -3885.103| 20| 21.0772|           0.3822| 1899.061|     0.0610|    0.1246|
|  3| -3877.990| 28| 24.2548|           0.0914| 1874.409|     0.0595|    0.1250|
|  4| -3873.144| 35| 22.5673|           0.0685| 1840.482|     0.0533|    0.1279|

Computational issues seem to arise as $K$ increases.
The best choice of $K$ according to the DIC criterion is $K=2$. We run again a longer MCMC chain fixing $K=2$. In particular, the number of simulations is set to 1000000, the thinning is set to 5 and the burn-in period is of 250000 observations. As an example, we report the Monte Carlo markov Chains for 3 parameters: $a_{11}, a_{a21}, a_{22}$ and $\sigma$.



```r
a           <- as.mcmc(   cbind(as.matrix(r2$A[,1:3,1]),r2$sigma2  )  )
colnames(a) <- c(rep(paste0("a[",c(11,21,31),"]")),"sigma")
aa          <- ggs(a)
a1          <- ggs_histogram(aa)+theme_bw()
a2          <- ggs_running(aa)+theme_bw() + theme(legend.position="none")
grid.arrange(a1, a2,  ncol=3, nrow=1)
```

![](tenStock_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

The prediction on the test dataset is performed. We then have a look at the resulting plots for a qualitative assessment. The model seems to capture the mean of the process.


```r
data.plot_DAT                   <- melt(DAT); colnames(data.plot_DAT) <- c("Stock","Time","Value")
data.plot_DAT$In_sample         <- factor(rep(1:n_t %in% time_grid,n_l)); 
levels(data.plot_DAT$In_sample) <- c("No","Yes")
data.plot_DAT$Pred_test         <- c(predict(r2, DAT_train, 1:n_t))

colnames(data.plot_DAT)[5]      <- paste("K=",2)
data.plot_DAT2                  <- melt(data.plot_DAT,id.vars = c("Stock","Time","Value","In_sample"))

ggplot(data=data.plot_DAT2 , aes(x=Time,y=Value)) + geom_point(alpha=0.25) + facet_wrap(   ~  Stock,ncol=2 ) + geom_line(aes(y=value),col=4, size=1.2) + ylab("Price")+ ggtitle("Ten Stocks")+ theme_bw()
```

![](tenStock_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Let us have a look to the raster plot which describes the estimated covariance matrix. Notice that this estimation is built averaging each component of the matrix along the simulated chain.




