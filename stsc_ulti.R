#stsc_ulti.py

compute_affinity_mat <- function(dist_mat, sigma = 1, local_scale = FALSE, k = 7) {
    sigma_sqr <- sigma^2
    
    if (local_scale) {
        kth_neighbor <- apply(dist_mat, 2, function(x) sort(x)[k])
        sigma_sqr <- tcrossprod(kth_neighbor)
    }
    
    # computing the weighted adjacency matrix
    W <- -dist_mat^2 / sigma_sqr
    W[is.nan(W)] <- 0
    W <- exp(W)
    diag(W) <- 0
    W
}

get_min_max <- function(values , min_n_cluster , max_n_cluster){
    args <- list(min_n_cluster, max_n_cluster)
    if (is.null(args[['min_n_cluster']])) {
        min_n_cluster <- 2
    }
    if (is.null(args[['max_n_cluster']])) {
        max_n_cluster <- sum(values[which(values >0)])
    }
    if (max_n_cluster < 2){
        max_n_cluster <- 2
    }
    if (min_n_cluster > max_n_cluster ){
        stop("min_n_cluster should be smaller than max_n_cluster")
    }
    list(min_n_cluster , max_n_cluster)
}


affinity_to_eigen2 <- function(W) {
  eps <- .Machine$double.eps
  degs <- colSums(W)
  D <- diag(degs)
  degs[degs == 0] <- eps
  diag(D) <- 1 / sqrt(degs)
  L <- D %*% W %*% D
  eigen(L)
}



refomat_result <- function(cluster_label , n){
    zip_data <- data.frame(cluster_label = cluster_label , number=1:n )
    zip_data[order(zip_data$cluster_label),]
    repeated <- as.numeric(names(table(cluster_label)[table(cluster_label) > 1]))
    grouped_idx = list()
    for ( i in 1:length(repeated)){
        grouped_idx[[i]] <- which(zip_data$cluster_label== repeated[i])
    }
    grouped_idx
}

