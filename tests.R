pth_filled_true <- "data/filled_true.csv"

true_wide <- read_csv(pth_filled_true, col_names = make_cl_nms(pth_filled_true), skip = 1) %>% 
  mutate(id = factor(id, levels = id))


true_long <- true_wide %>% 
  gather("occurrence", value = "time_raw", -id, na.rm = FALSE)

test_data <- true_long %>% 
  filter(!is.na(time_raw)) %>% 
  rename(time_raw_true = time_raw) %>% 
  inner_join(fitted_long, by = c("id", "occurrence")) %>%
  left_join(valid_filled_long, by = c("id", "occurrence")) %>% 
  filter(!is.na(time_centered)) %>% 
  rename(time_raw_hat = time_raw) %>% 
  left_join(main_long, by = c("id", "occurrence")) %>%
  mutate(time_raw = ifelse(!is.na(time_raw), time_raw, time_raw_hat))

Metrics::mdae(test_data$time_raw, test_data$time_raw_true)
