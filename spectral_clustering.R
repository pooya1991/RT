library(tidyverse)
library(RSpectra)
library(ClusterR)
library(speccalt)
# utils -------------------------------------------------------------------

silhouette_score <- function(clusters, dist) {
  sil <- cluster::silhouette(clusters, dist)
  mean(sil[, 3])
}


externalValidation <- function(y_true , y_pred){
  
  v1 <- external_validation(y_true , y_pred, method = "rand_index" , summary_stats=F)
  
  return(v1)
}


# analysis ----------------------------------------------------------------

X <- scan("data/X.csv", sep = ",")
X <- matrix(X, ncol = 2, byrow = TRUE)
Y <- scan("data/y_true.csv", sep = ",")
Y <- matrix(Y, ncol = 1)
Y <- Y+1

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
externalValidation(Y , cluster_predict)


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
cluster_predict = anocva::spectralClustering(W, 7)

# substituting anocva with speccalt package
kern <- local.rbfdot(X)
clust <- speccalt(kern, 7)

g1 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(clusters_true)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Ground truth simulated data : 7 clusters")

g2 <- ggplot(df, aes(x1, x2)) +
  geom_point(aes(color = factor(clust)), show.legend = FALSE, size = 3) +
  scale_color_viridis_d() +
  labs(x = NULL, y = NULL, title = "Spectral clustering result")

grid.arrange(g2, g1, ncol=2)

# cheking result by elbowmwthod , silhouette method============================================
## Clustering silhouette method
library(cluster)
pam_k7 <- pam(X , k=7)
pam_k7$silinfo$widths

sil_plot <- silhouette(pam_k7)
#plot(sil_plot)
pam_k7$silinfo$avg.width

# for different number of k
sil_width <- map_dbl(2:10 , function(k){
  mod <- pam( x=X , k=k)
  mod$silinfo$avg.width
})

sil_df <- data.frame(k= 2:10 , sil_width = sil_width)
#print(sil_df)

ggplot(sil_df , aes(x=k , y= sil_width)) + geom_line()+
  scale_x_continuous(breaks = 2:10)+ labs(title = "Silhouette score values vs Numbers of Clusters")

## Generating differnt kmean with different k - elbow method============================================
tot_within <- map_dbl(1:10 , function(k){
  model <- kmeans ( x= X , centers = k)
  model$tot.withinss
  
})
elbow_df <- data.frame( k = 1:10 , tot_within =tot_within )

print(elbow_df)

ggplot(elbow_df , aes(x=k , y= tot_within , color = "chartreuse"))+
  geom_line(show.legend = F , linetype = "dashed")+ geom_point(color="blue" , size=1)+
  labs(title = "The Elbow Method showing the optimal k")

## kmeans cluatering
model_k7 <- kmeans(X , centers = 7)
clust_k7 <- model_k7$cluster



