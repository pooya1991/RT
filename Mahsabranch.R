### Topology Overlap Matrix : In this method first i calculate the correlation between our variables.
# i use the correlation matrix to derive the adjacency matrix. there are two choices whether the 
# relation between Variables will be weighted(soft thresholding) or unweighted(hard thresholding).
# in the latter the resulting adjacency matrix will have just 0 , 1. 

install.packages("BiocManager")
BiocManager::install("WGCNA")
library(tidyverse)
library(viridis)

### the valid_centered_wide data has 3631 sample(IDs) and 99 featurs(occr)

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


# calculating adjacency matric for hard threshold (t0 = 0.7) by the choice of parameters

adj_mat1 <- simil_mat
size <- 3631
t0 <- 0.7

for (i in 1:size ){
    for (j in 1:size){
        if (adj_mat1[i,j] > t0 ) {
            adj_mat1[i,j] <- 1
        } else {
            adj_mat1[i,j] <- 0
        }
    }
}
diag(adj_mat1) <- 0

# calculating adjacency matrix with sigmoid function
alpha <- 10 
t0 <- 0.5
adj_mat2 <- 1/(1 + exp(-alpha*(simil_mat-t0)))
diag(adj_mat2) <- 0

# calculating adjacency matrix with power function(soft threshold), i checked for positive and symmetricity of
# similarity matrix and adjacency matix
# adjacency = power(similarity , \beta) ,\beta is the parameter should be chosen. 

A <- adjacency.fromSimilarity(simil_mat,
                         type = "unsigned",
                         power = 6)

k <- as.numeric(apply(A,2, sum))-1 
# standardize the connectivity 
Z.k <- scale (k)
max(Z.k)
min(Z.k)
thresholdZ.k <- -2.5
outliercolor <- ifelse(Z.k < thresholdZ.k, "black" , "red")
connectivity_color <- data.frame(numbers2colors(Z.k,signed = TRUE))
datcolor <- data.frame(outlier=outliercolor, connectivity_color)
IDsTree <- hclust(as.dist(1-A) , method = "average")
pdf("dendogram IDs.pdf")
plotDendroAndColors(IDsTree,colors = datcolor, groupLabels = names(datcolor), main ="Dendogram of IDs")
dev.off()



    
#Analysis of scale free topology for multiple soft thresholding powers
sft <- pickSoftThreshold.fromSimilarity(
    simil_mat,
    powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
    removeFirst = FALSE, nBreaks = 10, blockSize = 1000,
    moreNetworkConcepts=FALSE,
    verbose = 2, indent = 0)


newdata <- valid_centered_wide %>% 
    select(-id) %>%
    data.matrix %>%
    # transposing the matrix is needed because we want to compute cosine similarity between IDs
    # so we need IDs to be the columns of the matrix
    t() %>%
    `colnames<-`(valid_centered_wide$id)
pdf("boxplot of IDs.pdf")
boxplot(newdata , show.names=T , las= 2 ,  cex.axis = 0.5 , col=rainbow(10, alpha=0.5))
dev.off()