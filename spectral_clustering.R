library(tidyverse)
library(ggplot2)
library(RSpectra)

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


# Using function spectralclustering==================================
library(anocva)
library(gridExtra)
cluster_predict = spectralClustering(W, 7)

g1 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(clusters_true)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Ground truth simulated data : 7 clusters")

g2 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(cluster_predict)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Spectral clustering result")

grid.arrange(g2, g1, ncol=2)

# cheking result by elbowmwthod , hierarchiacal clustering 
## hierarchiacal clustering 
dist_data <- dist(X , method = 'euclidean')
hc_dist_data <- hclust(dist_data , method ='complete')

cluster <- cutree(hc_dist_data , h=5)
plot(hc_dist_data  , labels = F , hang = -1)
abline(h=5 , col = "red")

## Generating differnt kmean with different k - elbow method
library(purrr)
tot_within <- map_dbl(1:10 , function(k){
  model <- kmeans ( x= X , centers = k)
  model$tot.withinss
  
})
elbow_df <- data.frame( k = 1:10 , tot_within =tot_within )

print(elbow_df)

ggplot(elbow_df , aes(x=k , y= tot_within , color = "chartreuse"))+
  geom_line(show.legend = F , linetype = "dashed")+ geom_point(color="blue" , size=1)
## kmeans cluatering
model_k7 <- kmeans(X , centers = 7)
clust_k7 <- model_k7$cluster

Y <- Y+1
results <- as_tibble(cbind(cluster , clust_k7 , Y), .name_repair = "unique") %>% 
  set_names(c("hierarchiacal", "kmeans", "clusters_true"))

members <- rbind(hierarchiacal= table(cluster), kmeans=table(clust_k7) , cluster_true =table(Y) )

g3 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(cluster)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "hierarchiacal clustering result")

g4 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(clust_k7)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "kmeans clustering result")

grid.arrange(g1, g2, g3 , g4 , nrow=2 , ncol=2)
