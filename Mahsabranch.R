install.packages("BiocManager")
BiocManager::install("WGCNA")
library(tidyverse)
library(viridis)

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

#connectivity
k <- colSums(A) - 1
# standardize the connectivity , based on max , min values in connectivity, i define a threshold.
# values which are around minimum i called them outliers(black color) 
Z.k <- scale (k)
max(Z.k)
min(Z.k)
thresholdZ.k <- -2.5
outliercolor <- ifelse(Z.k < thresholdZ.k, "black" , "red")
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
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft <- pickSoftThreshold.fromSimilarity(
    simil_mat,
    powerVector = powers,
    removeFirst = FALSE, nBreaks = 10, blockSize = 1000,
    moreNetworkConcepts=FALSE,
    verbose = 2, indent = 0)

pdf("scale Independance.pdf")
par(mfrow=c(1,2))
cexl=0.5
plot(sft$fitIndices[,1] , -sign(sft$fitIndices[,3])*sft$fitIndices[,2] , xlab= "soft threshold" ,
     ylab = "scale free topology model fit" , type="n" , main ="scale Independence")
text(sft$fitIndices[,1] ,-sign(sft$fitIndices[,3])*sft$fitIndices[,2] , labels = powers , cex=0.5 ,
     col = "red")

plot(sft$fitIndices[,1] , sft$fitIndices[ , 5] , xlab = "soft threshold" , ylab = "mean connectivity",
     type = "n" , main = "mean connectivity")
text(sft$fitIndices[,1] , sft$fitIndices[,5] , labels = powers , cex=0.5 , col = "red")
dev.off()
