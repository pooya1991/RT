# template for checking the results of my code 
source("stsc_ulti.R")
source("stsc.R")
source("stsc_final.R")

x <- scan("base.csv", sep = ",") %>% matrix(ncol = 20)
distx <- as.matrix(dist(t(x)))

#affinity_matrix
W <- compute_affinity_mat(distx)
self_tuning_spectral_clustering(W)
