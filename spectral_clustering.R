library(tidyverse)
library(RSpectra)

# utils -------------------------------------------------------------------

silhouette_score <- function(clusters, dist) {
  sil <- cluster::silhouette(clusters, dist)
  mean(sil[, 3])
}

external_validation <- function(true_labels , clusters){
  tbl <- table(clusters, true_labels)
  
  tp_plus_fp <- sum(choose(rowSums(tbl), 2))
  tp_plus_fn <- sum(choose(colSums(tbl), 2))
  tp <- sum(choose(drop(tbl), 2))
  fp <- tp_plus_fp - tp
  fn <- tp_plus_fn - tp
  tn <- choose(sum(drop(tbl)), 2) - tp - fp - fn
  prod_comb <- (tp_plus_fp * tp_plus_fn) / choose(length(true_labels), 2)
  mean_comb <- (tp_plus_fp + tp_plus_fn) / 2
  adjusted_rand_index <- (tp - prod_comb) / (mean_comb - prod_comb)
  
  purity <- sum(apply(tbl, 1, max)) / length(true_labels)

  list(adjusted_rand_index = adjusted_rand_index, purity = purity)

}

compute_affinity_mat <- function(X , k){
  dist_mat <- dist(X) %>% as.matrix()
  kth_neighbor <- apply(dist_mat, 2, function(x) sort(x)[k])
  local_scale <- tcrossprod(kth_neighbor)
  
  # computing the weighted adjacency matrix
  W <- -dist_mat^2 / local_scale
  W[is.nan(W)] <- 0
  W <- exp(W)
  diag(W) <- 0
  W
}



eigen_value <- function(W , n_eigenval){
  deg <- colSums(W!= 0)
  D <- diag(1 / sqrt(deg))
  # Normalized laplacian matrix
  n = nrow(D)
  L <- diag(n) - D %*% W %*% D
  # eigenvalues with largest magnitude 
  eigenvalues <- eigen(N_L)$values
  eigenvectors <- eigen(N_L)$vectors
  plot(1:length(eigenvalues) ,  sort(eigenvalues , decreasing = F), main = "largest eigenvalues of matrix" , type="p" , col="blue" ,
       pch = 21 , bg="blue" )
  largest_gap <- abs(diff(eigenvalues))
  N <- n_eigenval
  ndx <- order(largest_gap, decreasing = T)[1:N]
  largest_gap[ndx]
  ndx
}


# analysis ----------------------------------------------------------------

X <- scan("data/X.csv", sep = ",")
X <- matrix(X, ncol = 2, byrow = TRUE)
Y <- scan("data/y_true.csv", sep = ",")
Y <- matrix(Y, ncol = 1)

# visualizing true clusters
df <- as_tibble(cbind(X, Y), .name_repair = "unique") %>% 
  set_names(c("x1", "x2", "clusters_true")) %>% 
  mutate(clusters_true = as.integer(clusters_true))

ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(clusters_true)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Ground truth simulated data : 7 clusters")


# spectral clustering using packages --------------------------------------

set.seed(942)
obj_clust <- kernlab::specc(X, centers = 7, kernel = "rbfdot", kpar = "automatic")

g1 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(clusters_true)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Ground truth clustering")

g2 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(obj_clust)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Spectral clustering results")

gridExtra::grid.arrange(g1 , g2, ncol = 2)

# internal validation
silhouette_score(obj_clust, dist(X))

# external validation
external_validation(drop(Y) , obj_clust)

# self tuned clustering ---------------------------------------------------


W <- compute_affinity_mat(X, 11)
LAV <- eigen_value(AF_MAT , n_eigenval = 5)
