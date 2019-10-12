library(tidyverse)
source("utils.R")

pth_filled_true <- "data/filled_true.csv"
pth_raw_data <- "data/for_Pooya_15.csv"
thresh_outlier <- 120
thresh_crop <- 60
n_clusts <- 10

true_wide <- read_csv(pth_filled_true, col_names = make_cl_nms(pth_filled_true), skip = 1) %>% 
    mutate(id = factor(id, levels = id))

true_long <- true_wide %>% 
    gather("occurrence", value = "time_raw", -id, na.rm = FALSE)

main_raw <- read_csv(pth_raw_data, col_names = make_cl_nms(pth_raw_data), skip = 1) %>% 
    mutate(id = factor(id, levels = id))

main_long <- main_raw %>% 
    gather("occurrence", value = "time_raw", -id, na.rm = FALSE)

main_long_aug <- main_long %>% 
    group_by(id) %>% 
    mutate(
        median_per_id = median(time_raw, na.rm = TRUE),
        time_centered = time_raw - median_per_id,
        time_raw      = ifelse(abs(time_centered) <= thresh_outlier, time_raw, NA_real_),
        time_centered = ifelse(abs(time_centered) <= thresh_outlier, time_centered, NA_real_),
        not_na_per_id = sum(!is.na(time_centered))
    ) %>% 
    ungroup()

valid_centered_long <- main_long_aug %>% 
    filter(not_na_per_id > 6) %>% 
    mutate(
        time_centered = ifelse(time_centered > thresh_crop, thresh_crop, time_centered),
        time_centered = ifelse(time_centered < -thresh_crop, -thresh_crop, time_centered)
    ) %>% 
    select(id, occurrence, time_centered)

valid_centered_wide <- valid_centered_long %>% 
    spread(occurrence, time_centered) %>% 
    arrange(id)

simil_mat <- valid_centered_wide %>% 
    select(-id) %>%
    data.matrix %>%
    t() %>%
    `colnames<-`(valid_centered_wide$id) %>%
    proxy::simil(method = "cosine", by_rows = FALSE) %>% 
    as.matrix()

simil_mat[near(abs(simil_mat), 1)] <- 0
diag(simil_mat) <- 1
simil_mat[is.na(simil_mat)] <- 0

W <- compute_affinity_mat(1 - simil_mat, local_scale = TRUE)
obj_eig <- affinity_to_eigen(W)
find_n_clusts(obj_eig)

eigen_vecs <- simil_to_eigen(simil_mat, n = 10, sigma = 3)
set.seed(1215)
res_clust <- e1071::cmeans(eigen_vecs, centers = n_clusts, m = 1.15)

ids_clusts <- tibble(
    id = valid_centered_wide$id,
    cluster = res_clust$cluster
)

fitter_factory <- fit_lasso_factory
fitter <- fitter_factory(valid_centered_wide, ids_clusts)
fitter_safe <- possibly(fitter, NA_real_)
future::plan(future::multiprocess)

fitted_long <- true_long %>% 
    filter(!is.na(time_raw)) %>% 
    select(id, occurrence) %>% 
    inner_join(filter(valid_centered_long, is.na(time_centered)), by = c("id", "occurrence")) %>% 
    mutate(time_centered = furrr::future_map2_dbl(id, occurrence, fitter_safe, .progress = TRUE))

valid_filled_long <- fitted_long %>% 
    left_join(select(main_long_aug, id, occurrence, median_per_id), by = c("id", "occurrence")) %>% 
    mutate(time_raw = time_centered + median_per_id) %>% 
    select(id, occurrence, time_raw)

test_data <- true_long %>% 
    rename(time_raw_true = time_raw) %>% 
    right_join(valid_filled_long, by = c("id", "occurrence")) %>% 
    rename(time_raw_hat = time_raw) %>% 
    left_join(main_long, by = c("id", "occurrence")) %>%
    mutate(time_raw = ifelse(!is.na(time_raw), time_raw, time_raw_hat))

Metrics::mdae(test_data$time_raw, test_data$time_raw_true)
Metrics::mae(test_data$time_raw, test_data$time_raw_true)

