# Model 1
Alessia Caponera, Francesco Denti, Tommaso Rigon, Andrea Sottosanti and Alan Gelfand  



First of all, we load the necessary data in memory. The data matrix `B` is composed by

- The observations of the first hemisphere (regions from 1 to 35).
- The `10`-th subject.
- First scan


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

The matrix `D` contains the euclidean distances among centroids. Again, we retained only the first `35` regions. In the `W` matrix we stored, instead the **proximity matrix** for the CAR prior.


```r
D  <- as.matrix(dist(subset(ROI, select=c(centroid_x,centroid_y,centroid_z))))[1:35,1:35]

index_D   <-  matrix(TRUE,nrow(D),nrow(D)); diag(index_D)<- FALSE
d_lower <- min(D[index_D])
d_upper <- max(D[index_D])

k1       <- 10^(d_lower/(d_upper - d_lower))
k2       <- log(k1)/d_lower
W       <- k1*exp(-k2*D); diag(W) <- 0
```


