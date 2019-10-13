# distance matrix 

data <-  valid_centered_wide %>% 
    select(-id) %>%
    data.matrix %>%
    # transposing the matrix is needed because we want to compute cosine similarity between IDs
    # so we need IDs to be the columns of the matrix
    t() %>%
    `colnames<-`(valid_centered_wide$id)


data_scaled <- scale(data)
distance <- dist(data_scaled, method = "manhattan")
distance <- as.matrix(distance)
W <- compute_affinity_mat(distance, sigma = 1, local_scale = TRUE, k = 7)
eig_obj <- affinity_to_eigen(W)
find_n_clusts(eig_obj , n=5)
