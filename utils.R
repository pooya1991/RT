# miscellaneous -----------------------------------------------------------

# this function calculates eigen vectors from a similarity matrix
simil_to_eigen <- function(simil_mat, n = 15, sigma = 3) {
    W <- exp(-(1 - simil_mat) / sigma)
    eps <- .Machine$double.eps
    degs <- rowSums(W)
    D <- diag(degs)
    L <- D - W
    degs[degs == 0] <- eps
    diag(D) <- 1 / (degs ^ 0.5)
    L <- D %*% L %*% D
    U <- eigen(L, symmetric = TRUE, only.values = FALSE)
    U <- U$vectors[, nrow(simil_mat):(nrow(simil_mat) - n + 1)]
    U <- diag(1 / sqrt(rowSums(U ^ 2))) %*% U
    U
}

# calc_stable_sd() calculates a stable standard deviation by trimming the observations in the upper and lower
# quartiles
calc_stable_sd <- function(x, na.rm = TRUE) {
    q <- quantile(x, c(.25, .75), na.rm = na.rm)
    st_sd <- (q[2] - q[1]) / (2 * .68)
    unname(st_sd)
}

make_cl_nms <- function(file_path) {
    cl_names <- readr::read_lines(file_path, n_max = 1) %>% 
        str_split(",") %>% unlist()
    cl_names[1] <- "id"
    cl_names[-1] <- paste0("occr", seq(100, length.out = length(cl_names[-1])))
    cl_names
}

# fitting functions -------------------------------------------------

# this function returns all IDs that are in the same cluster as the input ID
find_clust_ids <- function(id, clusts_lookup) {
    idx <- clusts_lookup == clusts_lookup[id]
    names(clusts_lookup)[idx]
}

# this function helps to throw out values of a vector that suspected as outliers
subset_b <- function(b, probs) {
    q <- quantile(b, probs = probs)
    b_trunc <- b[between(b, q[1], q[2])]
    sd_trunc <- sd(b_trunc)
    mean_trunc <- mean(b_trunc)
    thresh <- c(mean_trunc - 4 * sd_trunc, mean_trunc + 4 * sd_trunc)
    between(b, thresh[1], thresh[2])
}

# this is a function factory that generates functions used for fitting lasso models. The method use to implement
# this function is rather complicated. One has to have knowlege about function factories in order to understand
# how this function works
fit_lasso_factory <- function(wide_df, ids_and_clusts) {
    clusts_lookup <- ids_and_clusts$cluster
    names(clusts_lookup) <- as.character(ids_and_clusts$id)
    data_mat <- wide_df %>% 
        select(-id) %>% 
        data.matrix() %>% 
        `rownames<-`(wide_df$id)
    
    occrs <- colnames(data_mat)
    ids <- rownames(data_mat)
    data_lookup <- !is.na(data_mat)

    filled_transposed <- wide_df %>% 
        left_join(ids_and_clusts, by = "id") %>% 
        group_by(cluster) %>% 
        mutate_at(vars(-id, -cluster),
                  ~ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)) %>% 
        ungroup(cluster) %>% 
        select(-id, -cluster) %>% 
        data.matrix() %>%
        `rownames<-`(wide_df$id) %>% 
        t()
    rm(data_mat, ids_and_clusts, wide_df)
    
    function(id, occr) {
        # set.seed(as.integer(id))
        id <- as.character(id)
        sub_occrs <- union(occr, occrs[data_lookup[id,]])
        sub_ids <- union(
            id, 
            intersect(find_clust_ids(id, clusts_lookup = clusts_lookup),
                      ids[data_lookup[, occr]])
        )
        wt <- colSums(data_lookup[sub_ids, sub_occrs])[-1] + .01
        lasso_mat <- filled_transposed[sub_occrs, sub_ids]
        train_x <- lasso_mat[-1, -1]
        train_y <- lasso_mat[-1, 1]
        glmnet_cv_args <- list(
            x = train_x, y = train_y,
            nlambda = 50, grouped = FALSE, nfold = 5, type.measure = "mae",
            alpha = 1, thresh = 1e-03, weights = wt
        )
        glmnet_cv_obj <- tryCatch(do.call(glmnet::cv.glmnet, glmnet_cv_args), error = function(e) NULL)
        if (is.null(glmnet_cv_obj)) {
            best_lam <- 100
        } else {
            best_lam <- glmnet_cv_obj$lambda.1se
        }
        glmnet_args <- list(
            x = train_x, y = train_y,
            lambda = best_lam,
            alpha = 1, thresh = 1e-07, weights = wt
        )
        glmnet_mdl <- do.call(glmnet::glmnet, glmnet_args)
        predict(glmnet_mdl, newx = lasso_mat[1, -1, drop = F], type = "response") %>% 
          as.vector()
    }
}

# filling functions -------------------------------------------------------

fill_using_clust_factory <- function(clust_centroids_df) {
    centroids_mat <- clust_centroids_df %>%
        select(-cluster) %>% 
        data.matrix() %>% 
        t() %>% 
        `colnames<-`(clust_centroids_df$cluster)
    function(df) {
        vec <- as.vector(as.matrix(df))
        vec_centered <- vec - median(vec, na.rm = TRUE)
        which_clust <- apply(
            centroids_mat, 2,
            function(x) {coop::cosine(x, vec_centered, use = "complete.obs")}
        ) %>% abs() %>%  which.min() %>% names()
        chosen_centroid <- centroids_mat[, which_clust]
        dc <- median(vec, na.rm = TRUE)
        vec_filled <- ifelse(!is.na(vec), vec, dc + chosen_centroid)
        matrix(vec_filled, nrow = 1) %>% 
            as.tibble() %>% 
            `colnames<-`(colnames(df))
    }
}

fill_using_data_factory <- function(valid_centered_filled_df) {
    filled_wide <- valid_centered_filled_df %>% 
        spread(occurrence, time_centered) %>% 
        select(-id)
    centroid_vec <- filled_wide %>% 
        summarise_all(median) %>% 
        as.matrix() %>% as.vector() %>% 
        setNames(colnames(filled_wide))
    rm(valid_centered_filled_df, filled_wide)
    function(df) {
        vec <- as.vector(as.matrix(df))
        dc <- median(vec - centroid_vec, na.rm = TRUE)
        vec_filled <- dc + centroid_vec
        matrix(vec_filled, nrow = 1) %>% 
            as.tibble() %>% 
            `colnames<-`(colnames(df))
    }
}
