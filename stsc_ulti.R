#stsc_ulti.py


get_max_min <- function(values , min_n_cluster , max_n_cluster){
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


affinity_to_eigen <- function(W) {
    eps <- .Machine$double.eps
    degs <- colSums(W)
    D <- diag(degs)
    L <- D - W
    degs[degs == 0] <- eps
    diag(D) <- 1 / sqrt(degs)
    L <- D %*% L %*% D
    eigen(L)
}

# eigenvalues are sorted in ascending order
values <- srot(eig_obj$values , decreasing = FALSE)
vectors <- eig_obj$vectors

# n= length(cluster_label)
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

