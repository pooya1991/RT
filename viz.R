graphs_clusters <- vector("list", 10)  
for (i in 1:n_clusts) {
    graphs_clusters[[i]] <- valid_centered_long %>%
        filter(!is.na(time_centered)) %>% 
        left_join(ids_clusts1, by = "id") %>% 
        filter(cluster == i) %>% 
        ggplot(aes(factor(occurrence), time_centered)) +
        geom_line(aes(group = id, color = id)) +
        scale_colour_viridis_d(option = "inferno") +
        theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
        xlab("Occurrence") +
        ylab("Time") +
        ggtitle(paste0("Cluster ", i))
}