# install.packages(c("tidyverse", "e1071", "proxy", "future", "furrr", "glmnet", "coop"))
library(tidyverse)
source("utils.R")

# set parameters ----------------------------------------------------------

pth_raw_data <- "data/for_Pooya_15.csv"
thresh_outlier <- 120
thresh_crop <- 60
n_clusters <- 10

# import data -----------------------------------------------------------------------

# main_raw is a data-frame consists of the raw data without any modification. The first 
# column of main_raw contains the IDs and the rest of the columns are measurements for
# each occurrence. So the number of rows equals the number of IDs and the number of columns
# is equal to the number of occurrences plus one.
main_raw <- read_csv(pth_raw_data, col_names = make_cl_nms(pth_raw_data), skip = 1) %>% 
    mutate(id = factor(id, levels = id))

# main_long contains the same information as main_raw. But instead of a wide data-frame 
# we have a long data-frame
main_long <- main_raw %>% 
    gather("occurrence", value = "time_raw", -id, na.rm = FALSE)

# from main_long, we compute new variables and then augment main_long with these new variables
main_long_aug <- main_long %>% 
    group_by(id) %>% 
    mutate(
        # median of time per ID
        median_per_id = median(time_raw, na.rm = TRUE),
        # centering time by subtracting medians
        time_centered = time_raw - median_per_id,
        # values of time for which the centered version devtiates more than 120 units from 
        # the median are considered outliers
        time_raw      = ifelse(abs(time_centered) <= thresh_outlier, time_raw, NA_real_),
        time_centered = ifelse(abs(time_centered) <= thresh_outlier, time_centered, NA_real_),
        # the number of non-missing measurements for each ID. Note that the measurements
        # that we considered as outliers and substituted with NAs are considerered 
        # missing measurements 
        not_na_per_id = sum(!is.na(time_centered))
    ) %>% 
    ungroup()

# partition data ----------------------------------------------------------

# in this step we partition the data into 4 parts. For each part we use a different rule
# to fill out missing values

# excluded_2_6 consists of IDs for which the numbers of valid measurement are between 2 and 6
excluded_2_6 <- main_long_aug %>% 
    filter(between(not_na_per_id, 2, 6)) %>% 
    select(id, occurrence, time_raw) %>% 
    spread(occurrence, time_raw) %>% 
    arrange(id)

# excluded_1 consists of IDs for which there's only one valid measurement
excluded_1 <- main_long_aug %>% 
    filter(not_na_per_id == 1) %>% 
    select(id, occurrence, time_raw) %>% 
    spread(occurrence, time_raw) %>% 
    arrange(id)

# excluded_0 consists of IDs for which there's no valid measurement
excluded_0 <- main_long_aug %>% 
    filter(not_na_per_id == 0) %>% 
    select(id, occurrence, time_raw) %>% 
    spread(occurrence, time_raw) %>% 
    arrange(id)

# valid_centered_long consists of IDs for which we have more than 6 valid measurements.
# These are the IDs that we consider valid and we use them to build models.
valid_centered_long <- main_long_aug %>% 
    filter(
        not_na_per_id > 6
    ) %>% 
    mutate(
        # values of time_centered which deviate more than 60 units from the median of 
        # the corresponding ID are fixed at 60 (or -60)
        time_centered = ifelse(time_centered > thresh_crop, thresh_crop, time_centered),
        time_centered = ifelse(time_centered < -thresh_crop, -thresh_crop, time_centered)
    ) %>% 
    select(id, occurrence, time_centered)

# valid_centered_wide contains the same information as valid_centered_long but is in 
# a wide format
valid_centered_wide <- valid_centered_long %>% 
    spread(occurrence, time_centered) %>% 
    arrange(id)

# clustering --------------------------------------------------------------

# In this step we perform clustering. To do so, first, we calculate the similarity matrix
# from valid_centered_wide. We use valid_centered_wide because it contains the information
# we need and converting it to a matrix is straightforward
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

# we extract the first 10 eigen vectors that are returned by simil_to_eigen(). Note that these
# are not the eigen vectors corresponding to the biggest eigen values.
eigen_vecs <- simil_to_eigen(simil_mat, n = 10, sigma = 1.5)
set.seed(1215)
# the parameter centers defines the number of clusters and m defines the degree of fuzziness
res_clust <- e1071::cmeans(eigen_vecs, centers = n_clusters, m = 1.3)

# for each clustering we perform, we build a ids_clust data-frame. These data-frames consists
# of two variables. One is IDs and the other is the cluster to which each ID is assigned
ids_clusts <- tibble(
    id = valid_centered_wide$id,
    cluster = res_clust$cluster
)

# model fitting -----------------------------------------------------------

fitter_factory <- fit_lasso_factory
# fitter_factory() takes a data-frame and a clustering lookup-table and returns
# a function that fits lasso models and is tailored for that particular data-frame
fitter <- fitter_factory(valid_centered_wide, ids_clusts)
# a safe version of fitter() which returns NA in case of error
fitter_safe <- possibly(fitter, NA_real_)

# multiprocess is a kind of parallelization which uses multicore evaluation if supported (Linux)
# otherwise uses multisession evaluation (Windows). Note that multicore evaluation is faster
future::plan(future::multiprocess)
# fitted_long consists of the combination of IDs and occurrences in valid_Centered_long
# for which the measurements were missing and we fill time_centered column with an
# appropriate estimation
fitted_long <- valid_centered_long %>%
    filter(is.na(time_centered)) %>%
    # slice(1:20) %>%
    # future_map2_dbl() is the parallelized version of map2_dbl which in turn is a variant
    # of map() family in purrr package. To put it simple, for each ID and occurrence, 
    # we call the tailored fitter() that returns an estimation of centered time.
    mutate(time_centered = furrr::future_map2_dbl(id, occurrence, fitter_safe, .progress = TRUE))

# filling -----------------------------------------------------------------

# In this step we fill the 4 incomplete data-frames, using a different method for each one

# firstly, valid_centered_long is filled with the estimations obtained from the previous step
valid_centered_filled <- fitted_long %>% 
    bind_rows(filter(valid_centered_long, !is.na(time_centered))) %>% 
    left_join(ids_clusts, by = "id") %>% 
    group_by(cluster, occurrence) %>% 
    mutate(
        # after imputing the missing values in the last step, we may still have some NAs that corresponds to
        # observations for which fit_lasso() have failed to generate and estimation. These missing values
        # are imputed using the within clusters medians
        time_centered = ifelse(
            is.na(time_centered), median(time_centered, na.rm = TRUE), time_centered
        )
    ) %>% 
    ungroup(cluster, occurrence) %>% 
    select(id, occurrence, time_centered)

# valid_filled consists of IDs that we assumed to be valid. By adding centered version of time
# and the median that we computed for each id (and is included in main_long_aug), we obtain
# the non-centered verion of time (time_raw).
valid_filled_long <- valid_centered_filled %>% 
    left_join(select(main_long_aug, id, occurrence, median_per_id), by = c("id", "occurrence")) %>% 
    mutate(time_raw = time_centered + median_per_id) %>% 
    select(id, occurrence, time_raw)

valid_filled_wide <- valid_filled_long %>% 
    spread(occurrence, time_raw)

# each row of cluster-centroids corresponds to the center of a particular cluster
cluster_centriods <- valid_centered_filled %>% 
    left_join(ids_clusts, by = "id") %>% 
    spread(occurrence, time_centered) %>% 
    select(-id) %>% 
    group_by(cluster) %>% 
    summarise_all(mean) %>% 
    arrange(cluster)

# fill_using_clust() are functions returned from another function, fill_using_clust_factory().
# The latter function takes cluster_centroids as input and outputs a function that estimates
# missing values for a particular ID based on the cluster to which that ID is closest to.
# We use this method for IDs contained in excluded_2_6
fill_using_clust <- fill_using_clust_factory(cluster_centriods)

excluded_2_6_filled <- excluded_2_6 %>% 
    group_by(id) %>% 
    nest() %>% 
    mutate(data = map(data, ~fill_using_clust(.x))) %>% 
    unnest()

# fill_using_data() are functions returned form another function, fill_using_data_factory().
# The latter function takes filled verion of valid_centered_long and returns a function that
# estimates missing values for a particular ID based on the distance of the single measurement
# that ID has to the corresponding measurement in valid_Centered_filled. We use this method
# for IDs contained in excluded_1
fill_using_data <- fill_using_data_factory(valid_centered_filled)

excluded_1_filled <- excluded_1 %>% 
    group_by(id) %>%
    nest() %>% 
    mutate(data = map(data, ~fill_using_data(.x))) %>% 
    unnest()

# for IDs that have no valid measurement, we substitute data medians computed from valid_filled
data_medians <- valid_filled_wide %>% 
    select(-id) %>% 
    summarise_all(median)

excluded_0_filled <- excluded_0 %>% 
    select(id) %>%
    mutate(data = list(data_medians)) %>% 
    unnest()

# output construction -----------------------------------------------------

# In this step we combine all the 4 partitions together and construct the final output
wide_semi_final <- bind_rows(
    valid_filled_wide,
    excluded_2_6_filled,
    excluded_1_filled,
    excluded_0_filled
) %>% 
    arrange(id)

wide_final <- wide_semi_final %>% 
    gather("occurrence", "time_new", -id, na.rm = FALSE) %>% 
    # here we substitute the measurement that we considered as outliers initially or fixed at
    # tresh_crop with their original value
    left_join(main_long, by = c("id", "occurrence")) %>%
    mutate(time_raw = ifelse(!is.na(time_raw), time_raw, time_new)) %>% 
    select(id, occurrence, time_raw) %>% 
    spread(occurrence, time_raw)
