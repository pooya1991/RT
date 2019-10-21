# stsc final 
source("stsc_ulti.R")
source("stsc.R")       
       
self_tuning_spectral_clustering <- function(affinity, get_rotation_matrix, min_n_cluster=NULL, max_n_cluster=NULL){
    eig_obj <- affinity_to_eigen2(W)
    # eigenvalues by default are sorted in decreasing order 
    values <- eig_obj$values
    # vectors are correspodence to eigenvalues
    vectors <- eig_obj$vectors
    min_n_cluster <- get_min_max(values, min_n_cluster, max_n_cluster)[[1]]
    max_n_cluster <- get_min_max(values, min_n_cluster, max_n_cluster)[[2]]
    re <- list()
    k <- 0
    for ( c in min_n_cluster:max_n_cluster){
        # eigenvectors correspodence to largest eigenvalues
        x <- vectors[ , 1:c]
        cost <-  get_rotation_matrix(x, c)[[1]]
        r <-  get_rotation_matrix(x,c)[[2]]
        k <- k+1
        re [[k]] <- list(cost , x%*%r)
        cat(" number of cluster :", c , "\t", "cost :" , cost)
       
    COST <- re[order(sapply(re,'[[',1))][[1]][[1]]
    Z <- re[order(sapply(re,'[[',1))][[1]][[2]]
    return(reformat_result(apply( Z, 1, which.max), dim(Z)[1]))
    }}
     
self_tuning_spectral_clustering_np <- function(affinity, min_n_cluster=NULL, max_n_cluster=NULL){
    self_tuning_spectral_clustering(affinity, get_rotation_matrix, min_n_cluster, max_n_cluster)
}
        
        