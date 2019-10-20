# stsc final 
source(c("stsc_ulti.R" , "stsc.R"))
       
       
self_tuning_spectral_clustering <- function(affinity, get_rotation_matrix, min_n_cluster=NULL, max_n_cluster=NULL){
    eig_obj <- affinity_to_eigen(W)
    # eigenvalues by default are sorted in decreasing order 
    values <- eig_obj$values
    # vectors are correspodence to eigenvalues
    vectors <- eig_obj$vectors
    min_n_cluster <- get_min_max(values, min_n_cluster, max_n_cluster)[1]
    max_n_cluster <- get_min_max(values, min_n_cluster, max_n_cluster)[2]
    re <- list()
    k <- 0
    for ( c in min_n_cluster:max_n_cluster){
        # eigenvectors correspodence to largest eigenvalues
        x <- vectors[ , 1:c]
        cost <-  get_rotation_matrix(x, c)[1]
        r <-  get_rotation_matrix(x,c)[2]
        k <- k+1
        re [k] <- c(cost , x%*%r)
        cat(" number of cluster :", c , "\t", "cost :" , cost)
        
        
    }}
     
self_tuning_spectral_clustering_np <- function(affinity, min_n_cluster=NULL, max_n_cluster=NULL){
    self_tuning_spectral_clustering(affinity, get_rotation_matrix_np, min_n_cluster, max_n_cluster)
}
        
        