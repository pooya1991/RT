# template for checking the results of my code 
source("stsc_ulti.R")
source("stsc.R")
source("stsc_final.R")


base <- matrix(rnorm(100*4), ncol=4)
x <- list()
for (i in 1:5){
    x[[i]] <-  base * rnorm(1)
}
x <- do.call(cbind, x)
x <- x+matrix(rnorm(100*20), ncol=20)*0.01

affinity_matrix <- compute_affinity_mat (x)
self_tuning_spectral_clustering(affinity_matrix)