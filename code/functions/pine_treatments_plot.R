## Purpose: plot treatment prevalence vs. pinyon & ponderosa ba
## Project: PIJA_Management

pine_treatments_plot <- function()
{
  library(tidyverse)
  library(arrow)
  
  ## bring in data ready for modeling
  d <- read_parquet("data/model_data.parquet") %>% 
    select(pinyon_ba, pin_ba_log, ponderosa_ba, pipo_ba_log,
           prescribed_fire, mastication, thinning, harvest) %>% 
    ## pivot longer by treatment category
    pivot_longer(
      cols = c(prescribed_fire, mastication, thinning, harvest),
      names_to = "treatment",
      values_to = "treatment_prop"
    ) %>% 
    filter(treatment_prop > 0) %>% 
    ## pivot longer by pine species for plotting
    pivot_longer(
      cols = c(pinyon_ba, ponderosa_ba),
      names_to = "species",
      values_to = "basal_area"
    ) %>% 
    mutate(
      ## convert proportion to percent
      treatment_pct = treatment_prop * 100,
      ## clean up species labels
      species = case_when(
        species == "pinyon_ba" ~ "Pinyon",
        species == "ponderosa_ba" ~ "Ponderosa"
      )
    ) 
  
  ## plot treatment percent vs. basal area, colored by species
  treatment_levels <- c("prescribed_fire", "mastication", "thinning", "harvest")
  treatment_labels <- c("Prescribed Fire", "Mastication", "Thinning", "Harvest")
  
  p <- ggplot(d, aes(x = basal_area + 0.01, y = treatment_pct, color = species)) +
    geom_point(alpha = 0.1) +
    facet_grid(~ factor(treatment, levels = treatment_levels, labels = treatment_labels), scales = "free_y") +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.35, base = exp(1)),
                       labels = scales::label_number(accuracy = 1),
                       breaks = c(0, 1, 5, 15, 50)) +
    scale_color_manual(
      values = c("Pinyon" = "#A8B98F", 
                 "Ponderosa" = "#8B4513"),
      name = "Pine Species"
    ) +
    labs(
      x = "Basal Area (m²/ha)",
      y = "Landscape Percent Treated"
    ) +
    theme_bw() +
    theme(legend.position = "top") +
    guides(color = guide_legend(override.aes = list(alpha = 1)))
  
  ## add horizontal dashed lines at max treatment percent observed
  max_treatments <- d %>%
    group_by(treatment) %>%
    summarize(max_pct = max(treatment_pct, na.rm = TRUE))
  
  p <- p +
    geom_hline(data = max_treatments, aes(yintercept = max_pct),
               linetype = "dashed", 
               ) +
    geom_text(data = max_treatments, aes(x = 0.02, y = max_pct, 
                                         label = sprintf("Max: %.0f%%", max_pct)),
              vjust = -0.5, hjust = 0, 
              size = 3, inherit.aes = FALSE)
  
  ## save plot
  ggsave("Figures/Figure_2.png", plot = p, width = 10, height = 5)
  
}
