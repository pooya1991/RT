install.packages("BiocManager")
BiocManager::install("WGCNA")
library(tidyverse)
library(viridis)
library(WGCNA)

# Topology Overlap Matrix ------------------------------------- 
# In this method first we need to calculate the correlation between our variables, then we derive 
# similarity matrix from correlation. As we find the similarity matrix we will be capable of calculating
# adjacency matrix;in this step we will consider different thresholds. There are two choices whether the 
# relation between Variables will be weighted(soft thresholding) or unweighted(hard thresholding).
# in the latter case the resulting adjacency matrix will have just 0 , 1. 
# I am using the WGCNA package. The WGCNA as an analysis method is described in :
# Zhang B and Horvath S (2005) 
# "A General Framework for Weighted Gene Co-Expression Network Analysis, Statistical Applications 
# in Genetics and Molecular Biology": Vol. 4: No. 1, Article 17 PMID: 16646834.

# import data -----------------------------------------------
# the valid_centered_wide is a dataframe consists 3631 rows and 99 columns.The first 
# column contains the IDs and the rest of the columns are measurements for
# each occurrence. So the number of rows equals the number of IDs and the number of columns
# is equal to the number of occurrences plus one.

#Similarity Matrix

simil_mat <- valid_centered_wide %>% 
    select(-id) %>%
    data.matrix %>%
    # transposing the matrix is needed because we want to compute cosine similarity between IDs
    # so we need IDs to be the columns of the matrix
    t() %>%
    `colnames<-`(valid_centered_wide$id) %>%
    proxy::simil(method = "cosine", by_rows = FALSE) %>% 
    # the output of simil() is not a matrix. But turning it into a matrix is simple
    as.matrix()

# it's very unlikely that two ids have cosine similarity equal to 1 (or -1). But there are 
# many data points in the constructed similarity matrix which are equal to 1. This happens
# when two IDs have only two overlapping measurements. In this case we can't claim that these
# IDs are similar to each other. So whenever we see a perfect similarity, we substitute it
# with 0.
simil_mat[near(abs(simil_mat), 1)] <- 0
diag(simil_mat) <- 1
# the missing entries of the similarity matrix are due to IDs that have less than 2
# overlapping measurements. We assign 0 to these entries.
simil_mat[is.na(simil_mat)] <- 0


# calculating adjacency matrix with power function(soft threshold),I checked for positive and
# symmetricity of similarity and adjacency matrix


# calculating adjacency matrix with power function(soft threshold), i checked for positive and symmetricity of
# similarity matrix and adjacency matix
# adjacency = power(similarity , \beta) ,\beta is the parameter should be chosen. 

A <- adjacency.fromSimilarity(simil_mat,type = "unsigned",power = 6)

# this calculates the connectivity
k <- colSums(A) - 1
# standardize the connectivity. Based on max , min values in connectivity, i define a threshold.
Z.k <- scale (k)
max(Z.k)
min(Z.k)
# Designate outliers if their Z.k value is below the threshold
thresholdZ.k <- -2.5
#the color vector indicates outliers (red)
outliercolor <- ifelse(Z.k < thresholdZ.k, "red" , "black")
connectivity_color <- data.frame(numbers2colors(Z.k,signed = TRUE))
datcolor <- data.frame(outlier=outliercolor, connectivity_color)
IDsTree <- hclust(as.dist(1-A) , method = "complete")

ids_clusts <- cutree(IDsTree , h=0.8)
ids_clusts <- ids_clusts %>% 
    enframe("id", "cluster") %>% 
    mutate(id = as.integer(id)) %>% 
    arrange(id) %>% 
    mutate(id = factor(id))

# In order to get clear values at the bottom of dendrogram in clustering, I used the following way of 
# generating the dendrogram plot, but still it is messy !!!

dendoIDsTree <- as.dendrogram(IDsTree)

svg(width = 60)
par(mfrow=c(2,1))
plot(dendoIDsTree ,cex=0.3 , main = "Main , methoe = complete"  )
plot(cut(dendoIDsTree, h=0.8)$upper, cex = 0.3 , main = "upper tree with cut at h=0.8")
dev.off()


   
# pdf("dendogram IDs.pdf")
# plotDendroAndColors(IDsTree,colors = datcolor, groupLabels = c("outlier= black","connectivitycolor"),
#                     main ="Dendogram of IDs" , cex.dendroLabels=0.3 , cex.rowText = 0.2)
# dev.off() 


### a beautiful view of dendogram in circle layout---- it needs to be modified !!!!
# pdf("fanview.pdf")
# plot(as.phylo(IDsTree), type = "fan" , cex = 0.1 )
# dev.off()

    
# Analysis of scale free topology for multiple soft thresholding powers
# I am trying to find the optimum power which yeild a scale free of variables
# Choose a set of soft thresholding powers
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft <- pickSoftThreshold.fromSimilarity(
    simil_mat,
    powerVector = powers,
    removeFirst = FALSE, nBreaks = 10, blockSize = 1000,
    moreNetworkConcepts=FALSE,
    verbose = 2, indent = 0)

pdf("scale Independance.pdf")
#par(mfrow=c(1,2))
cexl=0.5
# SFT index as a function of different powers
plot(sft$fitIndices[,1] , -sign(sft$fitIndices[,3])*sft$fitIndices[,2] , xlab= "soft threshold(power)" ,
     ylab = "scale free topology model fit R^2" , type="n" , main ="scale Independence")
text(sft$fitIndices[,1] ,-sign(sft$fitIndices[,3])*sft$fitIndices[,2] , labels = powers , cex=0.5 ,
     col = "red")

# Mean connectivity as a function of different powers
# plot(sft$fitIndices[,1] , sft$fitIndices[ , 5] , xlab = "soft threshold" , ylab = "mean connectivity",
#      type = "n" , main = "mean connectivity")
# text(sft$fitIndices[,1] , sft$fitIndices[,5] , labels = powers , cex=0.5 , col = "red")
dev.off()
 
# SPECTRAL CLUSTERING  -----------------------------------------------------------
# similarity matrix and spectral clustering 
library(RSpectra)

# distance matrix 
Ids_dist <-  1 - simil_mat

## Based on the view of difference between eigenvalues----------------------------
# The optimal number of clusters by eigengap heuristic
nn <- 7
size <- dim(Ids_dist)[1]
# intialize the k-nearest neighbor matrix
knn_mat <- matrix(0, nrow = size , ncol = size)
neighbor_index <- matrix(0, nrow = size , ncol = nn)
# find the 7 nearest neighbors for each point, considering local scale
for (i in 1:size){
    neighbor_index[i,] <- order(Ids_dist[i,])[2:(nn + 1)]
    knn_mat[i,][neighbor_index] <- 1 
}

colnames(neighbor_index) <- seq(1 , 7 , 1)
row.names(neighbor_index) <- colnames(Ids_dist)
# k'th position is chosen( here I consider 1th nearest neighbor)
# I find the local_scale.The local_scale is the distance between each Id and its 1th nearest neighbor
local_scale1 <- c()
for (i in 1:size){
    local_scale1 [i] <- Ids_dist[i , neighbor_index[ i,1]]
}

local_scale1 <- matrix(local_scale1 , nrow = size , ncol = 1)

# calculating sigma i * sigma j
sigma <- local_scale1 %*% t(local_scale1)

# Affinity matrix by considering local_scale ( 1th nearest neighbor)
affinitymatirx <- exp(-(Ids_dist*Ids_dist)/sigma)
diag(affinitymatirx) <- 0
# A diagonal matrix which diagonal elements are degree of each variable in affinity matrix
deg <- colSums(affinitymatirx != 0)
D <- diag(1 / sqrt(deg))
# laplacian matrix
L <- diag(deg) - affinitymatirx
# Normalized laplacian matrix
N_L <- D %*% L %*% D

# calculating eigenvalues and eigenvectors of laplacian matrix
eigenvalues <- eigen(N_L)$values
largest_gap <- abs(diff(eigenvalues))


plot(1:10 , eigenvalues[1:10] , main = "largest eigenvalues of matrix" , type="p" , col="blue" ,
     pch = 21 , bg="blue" )


## Based on paper: self-tuning spectral clustering ---------------------------------------------------
nn <- 7
size <- dim(Ids_dist)[1]
# intialize the knn matrix
knn_mat <- matrix(0, nrow = size , ncol = size)
neighbor_index <- matrix(0, nrow = size , ncol = nn)
# find the 7 nearest neighbors for each point, considering local scale
for (i in 1:size){
    neighbor_index[i,] <- order(Ids_dist[i,])[2:(nn + 1)]
    knn_mat[i,][neighbor_index] <- 1 
}

colnames(neighbor_index) <- seq(1 , 7 , 1)
row.names(neighbor_index) <- colnames(Ids_dist)
# According to paper the 7th neighbor is chosen( here i consider 7th nearest neighbor)
# i find the local_scale.local_scale is the distance between each Id and its 7th nearest neighbor
local_scale <- c()
for (i in 1:size){
    local_scale [i] <- Ids_dist[i , neighbor_index[ i,nn]]
}


# deriving affinity matrix . I should divide the all columns of (distance)^2 by local_scale(vector)
affmat <- exp(-(Ids_dist*Ids_dist)/local_scale)
# degrees of vertices
g <- colSums( affmat !=0 )
 
# normalizatoin 
D_half <- diag(1 / sqrt(g) )

# normalized laplacian matrix
# subtracting normalized affinity matrix from Identity matrix 
normal_laplacian <- diag(size) - (D_half %*% affmat %*% D_half)

res = eigs_sym(normal_laplacian , k=15, which = "LM")  #"LM" is the default largest magnitude , 10 largest eigenvectors
eigenvalues <- res$values
eigenvectors <- res$vectors

# in this step people use kmeans with a MANUALLY chosen number of clusters, but I am following the
# algorithm of paper then there is another way to determine the optimum number of clusters
