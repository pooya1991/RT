library(tidyverse)

X <- scan("data/X.csv", sep = ",")
X <- matrix(X, ncol = 2, byrow = TRUE)
Y <- scan("data/y_true.csv", sep = ",")
Y <- matrix(Y, nrow = 1)

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
L <- D %*% W %*% D
