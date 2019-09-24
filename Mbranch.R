### Topology Overlap Matrix

### the valid_centered_wide data has 3631 sample(IDs) and 99 featurs(occr)

data =  valid_centered_wide
mat <- valid_centered_wide %>% 
    select(-id) %>%
    data.matrix %>%
    # transposing the matrix is needed because we want to compute cosine similarity between IDs
    # so we need IDs to be the columns of the matrix
    t() %>%
    `colnames<-`(valid_centered_wide$id) %>% as.matrix()
    

simil_mat <- (1+cor(mat))/2

# Clustering silhouette method
library(cluster)
pam_k15 <- pam(data , k=15)
pam_k15$silinfo$widths

sil_plot <- silhouette(pam_k15)

pdf("silhouette15.pdf")
plot(sil_plot)
dev.off()


pam_k10$silinfo$avg.width

ids_clustering <- tibble(
    id = valid_centered_wide$id,
    cluster = pam_k10$clustering
)

head(data)
