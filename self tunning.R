# distance matrix 
Ids_dist <-  1 - simil_mat
data <-  valid_centered_wide %>% 
    select(-id) %>%
    data.matrix %>%
    # transposing the matrix is needed because we want to compute cosine similarity between IDs
    # so we need IDs to be the columns of the matrix
    t() %>%
    `colnames<-`(valid_centered_wide$id)

A <- compute_affinity_mat(data , k=7)
degree <- colSums(A)
De <- diag(degree)
diag(De) <- 1 / (degree ^ 0.5)

# Normalize affinity matrix 
NL <-  De %*% A %*% De
eig_obj <- affinity_to_eigen(NL)
# i want to use similarity transformation to make NL matrix diagonal
LC <- eig_obj$vectors[ ,1:5]




