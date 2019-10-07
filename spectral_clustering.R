library(tidyverse)
library(RSpectra)
library(ClusterR)
library(speccalt)
library(caret)

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

# internal validation
silhouette_score(drop(Y), dist(X))


# externalvalidation
external_validation (Y , clust2)


k <- 11
dist_mat <- dist(X) %>% as.matrix()
kth_neighbor <- apply(dist_mat, 2, function(x) sort(x)[k])
local_scale <- tcrossprod(kth_neighbor)

# computing the weighted adjacency matrix
W <- -dist_mat^2 / local_scale
W[is.nan(W)] <- 0
W <- exp(W)
diag(W) <- 0

# computing the degree matrix
eps <- .Machine$double.eps
degs <- colSums(W)
degs[degs == 0] <- eps
D <- diag(1 / sqrt(degs))

# computing the normalized laplacian matrix
n = nrow(D)
L <- diag(n) - D %*% W %*% D

# eigenvalues with largest magnitude 
##res = eigs_sym(L , k=k, which = "LM") !!!!
eigenvalues <- eigen(L)$values
eigenvectors <- eigen(L)$vectors

pdf("largest eigenvalues of matrix.pdf")
plot( 1:length(eigenvalues) , sort(eigenvalues , decreasing = F) , main = "largest eigenvalues of matrix" , type="p",
      col="blue" ,pch = 16 , bg="red" )
dev.off()

# recommended optimal number of clusters 
index_largest_gap <- which.max(abs(diff(sort(eigenvalues , decreasing = F))))


# Using function spectralclustering ================================================
library(anocva)
library(gridExtra)
cluster_predict = anocva::spectralClustering(W, 7)

# substituting anocva with speccalt package
kern <- local.rbfdot(X)
clust <- speccalt(kern, 7)

# kernlab package
# class data is matrix
clust2 <- kernlab::specc(X, centers=7, kernel = "rbfdot", kpar = "automatic")[1:length(Y)]


g1 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(clusters_true)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Ground truth simulated data : 7 clusters")

g2 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(clust)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Spectral clustering result by specclat package")

g3 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(clust2)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Spectral clustering result by kernlab package ")



grid.arrange(g1 , g2 , g3, ncol=3)

# confusion matrix =======================================================================
library(caret)
reference <- factor(Y)
data <- factor(clust2)
confusion_Matrix <- caret::confusionMatrix(data, reference, positive = NULL,
                dnn = c("Prediction", "Reference"), prevalence = NULL,
                mode = "sens_spec")


