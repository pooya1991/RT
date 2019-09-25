### Topology Overlap Matrix : In this method first i calculate the correlation between our variables.
# i use the correlation matrix to derive the adjacency matrix. there are two choices whether the 
# relation between Variables will be weighted(soft thresholding) or unweighted(hard thresholding).
# in the latter the resulting adjacency matrix will have just 0 , 1. 


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
simil_mat <- (1+ correlation)/2


# all the elements are NA

# Clustering silhouette method
library(cluster)
pam_k15 <- pam(data , k=15)
pam_k10$silinfo$widths

sil_plot <- silhouette(pam_k15)

pdf("silhouette15.pdf")
plot(sil_plot)
dev.off()


pam_k10$silinfo$avg.width

