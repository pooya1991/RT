# template for checking the results of my code 
source("stsc_ulti.R")
source("stsc.R")
source("stsc_final.R")

W <- scan("similarity_matrix.csv", sep = ",") %>% matrix(ncol = 20)
self_tuning_spectral_clustering(W)
