### Topology Overlap Matrix : In this method first i calculate the correlation between our variables.
# i use the correlation matrix to derive the adjacency matrix. there are two choices whether the 
# relation between Variables will be weighted(soft thresholding) or unweighted(hard thresholding).
# in the latter the resulting adjacency matrix will have just 0 , 1. 

install.packages("BiocManager")
BiocManager::install("WGCNA")

### the valid_centered_wide data has 3631 sample(IDs) and 99 featurs(occr)

data =  valid_centered_wide
mat_data <- valid_centered_wide %>% 
    select(-id) %>%
    data.matrix %>%
    # transposing the matrix is needed because we want to compute cosine similarity between IDs
    # so we need IDs to be the columns of the matrix
    t() %>%
    `colnames<-`(valid_centered_wide$id) %>% as.matrix()

correlation <- cor(mat_data , use = "pairwise.complete.obs")
# this correlation still has NA it stems from some variables which has zero variance, for example if 
# get correlation cbind(a=runif(10),b=rep(1,10)) i will face the error" the standard deviation is zero"
# i want to remove zero variance predictor in df.

dfr <- mat_data[, ! apply(mat_data , 2 , function(x) sd(x, na.rm = TRUE)==0 ) ]

all.equal(dfr , mat_data)


simil_mat <- (1+ correlation)/2 # still NA !!!!!!!!!!!
diag(simil_mat) <- 0

any(is.na(simil_mat))
sum(is.na(simil_mat))
max(simil_mat , na.rm = T) 
min(simil_mat , na.rm = T)

# calculating adjacency matric for hard threshold by the choice of parameters

adj_mat1 <- si

# calculating adjacency matrix with sigmoid function 
alpha <- 10 
t0 <- 0.5
adj_mat <- 1/(1 + exp())

# DBSCAN
library(fpc)
library(ggplot2)
library(factoextra)
set.seed(123)
db <- fpc::dbscan(mat_data, eps = 0.15, MinPts = 5 )

fviz_cluster(db, data = df, stand = FALSE,
             ellipse = FALSE, show.clust.cent = FALSE,
             geom = "point",palette = "jco", ggtheme = theme_classic())

print(db)
print.dbscan()

# Determining the optimal eps value
par(mfrow = c(1,3))
dbscan::kNNdistplot(mat_data, k =  5)

# Clustering silhouette method
library(cluster)
pam_k15 <- pam(data , k=15)
pam_k10$silinfo$widths

sil_plot <- silhouette(pam_k15)

pdf("silhouette15.pdf")
plot(sil_plot)
dev.off()


pam_k10$silinfo$avg.width

