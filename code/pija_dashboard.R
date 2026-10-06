## Purpose: Pinyon Jay data management and analysis dashboard
## Project: PIJA_Management

## plot the prevalence of different treatment types with respect to pinyon and ponderosa pine density
source("code/functions/pine_treatments_plot.R")
pine_treatments_plot()

## run a brms model; takes about 10 hrs
source("code/functions/fit_zip_brms.R")
m <- fit_zip_brms(
  data_path = "data/intermediate/pinjay_pines_ba.parquet",
  output_path = "models/pinjay_zip.rds", 
  return_model = TRUE,
  use_approximation = FALSE,
  iter = 3000,
  chains = 6
)


## if not running, read in
fit <- readRDS("models/pinjay_zip.rds")

## pull out list of parameters (fixed effects only)
source('code/functions/list_available_effects.R')
list_available_effects(fit)

## extract and plot zip effects
source('code/functions/extract_and_plot_zip_effects.R')
extract_and_plot_zip_effects(
  model = fit,
  effects_to_include = rev(c(
    "year_scaled",
    "elevation_30m_median_scaled",
    # "Ielevation_30m_median_scaledE2",
    # "tcc2001_scaled",
    "propwater_scaled",
    # "day_of_year_scaled",
    # "Iday_of_year_scaledE2",
    # "cci_scaled",
    # "effort_hrs_scaled",
    # "num_observers_scaled",
    # "is_mobile",
    # "is_mobile:effort_distance_km_scaled",
    "pin_ba_log_scaled",
    "Ipin_ba_log_scaledE2",
    "pipo_ba_log_scaled",
    "Ipipo_ba_log_scaledE2",
    "prescribed_fire_bin:prescribed_fire_scaled",
    "prescribed_fire_bin:prescribed_fire_scaled:prescribed_fire_yrsince_scaled",
    "pin_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled",
    "pipo_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled",
    "mastication_bin:mastication_scaled",
    "mastication_bin:mastication_scaled:mastication_yrsince_scaled",
    "pin_ba_log_scaled:mastication_bin:mastication_scaled",
    "pipo_ba_log_scaled:mastication_bin:mastication_scaled",
    "thinning_bin:thinning_scaled",
    "thinning_bin:thinning_scaled:thinning_yrsince_scaled",
    "pin_ba_log_scaled:thinning_bin:thinning_scaled",
    "pipo_ba_log_scaled:thinning_bin:thinning_scaled",
    "harvest_bin:harvest_scaled",
    "harvest_bin:harvest_scaled:harvest_yrsince_scaled",
    "pin_ba_log_scaled:harvest_bin:harvest_scaled",
    "pipo_ba_log_scaled:harvest_bin:harvest_scaled"
    # "clearcut_bin:clearcut_scaled",
    # "clearcut_bin:clearcut_scaled:clearcut_yrsince_scaled",
    # "pin_ba_log_scaled:clearcut_bin:clearcut_scaled",
  )),
  effect_labels = c(
    "year_scaled" = "Year",
    "elevation_30m_median_scaled" = "Elevation",
    "Ielevation_30m_median_scaledE2" = "Elevation²",
    "pin_ba_log_scaled" = "Pinyon Pine BA",
    "Ipin_ba_log_scaledE2" = "Pinyon Pine BA²",
    "pipo_ba_log_scaled" = "Ponderosa Pine BA",
    "Ipipo_ba_log_scaledE2" = "Ponderosa Pine BA²",
    # "tcc2001_scaled" = "Tree Canopy Cover",
    "propwater_scaled" = "Proportion Water",
    "prescribed_fire_bin:prescribed_fire_scaled" = "Prescribed Fire",
    "prescribed_fire_bin:prescribed_fire_scaled:prescribed_fire_yrsince_scaled" = "Prescribed Fire x Time",
    "pin_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled" = "Prescribed Fire x Pinyon BA",
    "pipo_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled" = "Prescribed Fire x Ponderosa BA",
    "mastication_bin:mastication_scaled" = "Mastication",
    "mastication_bin:mastication_scaled:mastication_yrsince_scaled" = "Mastication x Time",
    "pin_ba_log_scaled:mastication_bin:mastication_scaled" = "Mastication x Pinyon BA",
    "pipo_ba_log_scaled:mastication_bin:mastication_scaled" = "Mastication x Ponderosa BA",
    "thinning_bin:thinning_scaled" = "Thinning",
    "thinning_bin:thinning_scaled:thinning_yrsince_scaled" = "Thinning x Time",
    "pin_ba_log_scaled:thinning_bin:thinning_scaled" = "Thinning x Pinyon BA",
    "pipo_ba_log_scaled:thinning_bin:thinning_scaled" = "Thinning x Ponderosa BA",
    "harvest_bin:harvest_scaled" = "Harvest",
    "harvest_bin:harvest_scaled:harvest_yrsince_scaled" = "Harvest x Time",
    "pin_ba_log_scaled:harvest_bin:harvest_scaled" = "Harvest x Pinyon BA",
    "pipo_ba_log_scaled:harvest_bin:harvest_scaled" = "Harvest x Ponderosa BA",
    # "clearcut_bin:clearcut_scaled" = "Clearcut",
    # "clearcut_bin:clearcut_scaled:clearcut_yrsince_scaled" = "Clearcut x Time",
    # "pin_ba_log_scaled:clearcut_bin:clearcut_scaled" = "Clearcut x Pinyon BA",
    "day_of_year_scaled" = "Day of Year",
    "Iday_of_year_scaledE2" = "Day of Year²",
    "cci_scaled" = "CCI",
    "effort_hrs_scaled" = "Effort Hours",
    "num_observers_scaled" = "Number of Observers",
    "is_mobile" = "Mobile Checklist",
    "is_mobile:effort_distance_km_scaled" = "Mobile Distance"
  ),
  show_component_prefix = F,
  include_zi_effects = F,
  csv_output = 'Data/Results/pinjay_zip_effects.csv',
  plot_output = 'Figures/pinjay_effects_main.png',
  # xlim = c(-1, .15),
  plot_height = 6,
  plot_width = 12
)

extract_and_plot_zip_effects(
  model = fit,
  effects_to_include = rev(c(
    "day_of_year_scaled",
    "Iday_of_year_scaledE2",
    "effort_hrs_scaled",
    "num_observers_scaled",
    "is_mobile",
    "cci_scaled",
    "is_mobile:effort_distance_km_scaled"
  )),
  effect_labels = c(
    "day_of_year_scaled" = "Day of Year",
    "Iday_of_year_scaledE2" = "Day of Year²",
    "effort_hrs_scaled" = "Effort Hours",
    "num_observers_scaled" = "Number of Observers",
    "is_mobile" = "Mobile Checklist",
    "is_mobile:effort_distance_km_scaled" = "Mobile Distance",
    "cci_scaled" = "CCI"
  ),
  show_component_prefix = T,
  include_zi_effects = T,
  csv_output = 'Data/Results/pinjay_zip_effects.csv',
  plot_output = 'Figures/pinjay_effects_si.png',
  # xlim = c(-1, .15),
  plot_height = 4,
  plot_width = 12
)

library(tidyverse)
fit <- readRDS("models/pinjay_zip.rds")
original_data <- arrow::read_parquet("data/intermediate/pinjay_pines_ba.parquet") %>% 
  mutate(is_mobile = if_else(is_stationary, 0L, 1L))

source("code/functions/marginal_effects_plots.R")

results <- plot_zip_marginal_effects(
  model = fit,
  original_data = original_data,
  predictors = c("year_scaled", "day_of_year_scaled", 
                 "effort_hrs_scaled",
                 "num_observers_scaled", "cci_scaled", 
                 "elevation_30m_median_scaled", 
                 "propwater_scaled",
                 "pin_ba_log_scaled",
                 "pipo_ba_log_scaled"),
  predictor_labels <- 
    c("Year", "Day of Year", "Effort (hrs)", "# Observers", 
      "CCI", "Elevation (m)", "Prop. Water", 
      "Pinyon BA (log)", "Ponderosa BA (log)"),
  n_draws = 500  # faster; use NULL for all draws
)

results$combined_plot
ggsave("figures/marginal_effects_linear.png", results$combined_plot, width = 10, height = 12)


results <- plot_zip_interaction_effects(
  model = fit,
  original_data = original_data,
  treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
  yrsince_levels = c(1, 5, 10),
  pinba_quantiles = c(0.2, 0.5, 0.8),
  pipoba_quantiles = c(0.2, 0.5, 0.8)
)

# Save the two figures separately
ggsave("figures/treatment_yrsince.png", results$combined_plot_yrsince, width = 14, height = 4)
ggsave("figures/treatment_ba.png", results$combined_plot_ba, width = 14, height = 8)

## flipping the interaction plots
## resource if needed
source("code/functions/marginal_effects_plots.R")
results_flipped <- plot_zip_interaction_effects_flipped(
  model = fit,
  original_data = original_data,
  treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
  prop_treated_levels = c(0, 0.1, 0.5, 0.8),  
  yrsince_range = c(1, 10),
  n_draws = 500
)


# Save the time-based figure
ggsave("figures/treatment_over_time.png", results_flipped$combined_plot_yrsince,
       width = 5, height = 8, dpi = 300)

# Save the BA-based figure (2 rows)
ggsave("figures/treatment_across_ba.png", results_flipped$combined_plot_ba,
       width = 12, height = 6, dpi = 300)


#### Seasonal model comparison ####
## let's also run models for the breeding and non-breeding seasons
## These take bout 20 hrs to run
source("code/functions/fit_zip_brms.R")
m_breeding <- fit_zip_brms(
  data_path = "data/intermediate/pinjay_pines_ba.parquet",
  output_path = "models/pinjay_zip_breeding.rds",
  return_model = TRUE,
  use_approximation = FALSE,
  season = "breeding",
  iter = 3000,
  chains = 6
)
source("code/functions/fit_zip_brms.R")
m_nonbreeding <- fit_zip_brms(
  data_path = "data/intermediate/pinjay_pines_ba.parquet",
  output_path = "models/pinjay_zip_nonbreeding.rds",
  return_model = TRUE,
  use_approximation = FALSE,
  season = "nonbreeding",
  iter = 3000,
  chains = 6
)

source("code/functions/compare_zip_effects_multimodel.R")

results <- compare_zip_effects_multimodel(
  models = list(
    # "Full Year" = readRDS("models/pinjay_zip.rds"),
    "Breeding" = readRDS("models/pinjay_zip_breeding.rds"),
    "Non-breeding" = readRDS("models/pinjay_zip_nonbreeding.rds")
  ),
  csv_output = "data/results/seasonal_comparison_main.csv",
  plot_output = "figures/seasonal_comparison_main.png",
  model_colors = c(
    "Breeding" = "#C4956A", 
    "Non-breeding" = "#404040"
  ),
  effect_groups = list(
    "Environmental" = rev(c("year_scaled",
                            "elevation_30m_median_scaled",
                            "propwater_scaled",
                            "pin_ba_log_scaled",
                            "Ipin_ba_log_scaledE2",
                            "pipo_ba_log_scaled",
                            "Ipipo_ba_log_scaledE2")),
    "Prescribed Fire" = rev(c(
      "prescribed_fire_bin:prescribed_fire_scaled",
      "prescribed_fire_bin:prescribed_fire_scaled:prescribed_fire_yrsince_scaled",
      "pin_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled",
      "pipo_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled"
    )),
    "Mastication" = rev(c(
      "mastication_bin:mastication_scaled",
      "mastication_bin:mastication_scaled:mastication_yrsince_scaled",
      "pin_ba_log_scaled:mastication_bin:mastication_scaled",
      "pipo_ba_log_scaled:mastication_bin:mastication_scaled"
    )),
    "Thinning" = rev(c(
      "thinning_bin:thinning_scaled",
      "thinning_bin:thinning_scaled:thinning_yrsince_scaled",
      "pin_ba_log_scaled:thinning_bin:thinning_scaled",
      "pipo_ba_log_scaled:thinning_bin:thinning_scaled"
    )),
    "Harvest" = rev(c(
      "harvest_bin:harvest_scaled",
      "harvest_bin:harvest_scaled:harvest_yrsince_scaled",
      "pin_ba_log_scaled:harvest_bin:harvest_scaled",
      "pipo_ba_log_scaled:harvest_bin:harvest_scaled"
    ))
  ),
  effects_to_include = rev(c(
    "year_scaled",
    "elevation_30m_median_scaled",
    "propwater_scaled",
    "pin_ba_log_scaled",
    "Ipin_ba_log_scaledE2",
    "pipo_ba_log_scaled",
    "Ipipo_ba_log_scaledE2",
    "prescribed_fire_bin:prescribed_fire_scaled",
    "prescribed_fire_bin:prescribed_fire_scaled:prescribed_fire_yrsince_scaled",
    "pin_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled",
    "pipo_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled",
    "mastication_bin:mastication_scaled",
    "mastication_bin:mastication_scaled:mastication_yrsince_scaled",
    "pin_ba_log_scaled:mastication_bin:mastication_scaled",
    "pipo_ba_log_scaled:mastication_bin:mastication_scaled",
    "thinning_bin:thinning_scaled",
    "thinning_bin:thinning_scaled:thinning_yrsince_scaled",
    "pin_ba_log_scaled:thinning_bin:thinning_scaled",
    "pipo_ba_log_scaled:thinning_bin:thinning_scaled",
    "harvest_bin:harvest_scaled",
    "harvest_bin:harvest_scaled:harvest_yrsince_scaled",
    "pin_ba_log_scaled:harvest_bin:harvest_scaled",
    "pipo_ba_log_scaled:harvest_bin:harvest_scaled"
  )),
  effect_labels = c(
    "year_scaled" = "Year",
    "elevation_30m_median_scaled" = "Elevation",
    "pin_ba_log_scaled" = "Pinyon Pine BA",
    "Ipin_ba_log_scaledE2" = "Pinyon Pine BA²",
    "pipo_ba_log_scaled" = "Ponderosa Pine BA",
    "Ipipo_ba_log_scaledE2" = "Ponderosa Pine BA²",
    "propwater_scaled" = "Proportion Water",
    "prescribed_fire_bin:prescribed_fire_scaled" = "Prescribed Fire",
    "prescribed_fire_bin:prescribed_fire_scaled:prescribed_fire_yrsince_scaled" = "Prescribed Fire x Time",
    "pin_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled" = "Prescribed Fire x Pinyon BA",
    "pipo_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled" = "Prescribed Fire x Ponderosa BA",
    "mastication_bin:mastication_scaled" = "Mastication",
    "mastication_bin:mastication_scaled:mastication_yrsince_scaled" = "Mastication x Time",
    "pin_ba_log_scaled:mastication_bin:mastication_scaled" = "Mastication x Pinyon BA",
    "pipo_ba_log_scaled:mastication_bin:mastication_scaled" = "Mastication x Ponderosa BA",
    "thinning_bin:thinning_scaled" = "Thinning",
    "thinning_bin:thinning_scaled:thinning_yrsince_scaled" = "Thinning x Time",
    "pin_ba_log_scaled:thinning_bin:thinning_scaled" = "Thinning x Pinyon BA",
    "pipo_ba_log_scaled:thinning_bin:thinning_scaled" = "Thinning x Ponderosa BA",
    "harvest_bin:harvest_scaled" = "Harvest",
    "harvest_bin:harvest_scaled:harvest_yrsince_scaled" = "Harvest x Time",
    "pin_ba_log_scaled:harvest_bin:harvest_scaled" = "Harvest x Pinyon BA",
    "pipo_ba_log_scaled:harvest_bin:harvest_scaled" = "Harvest x Ponderosa BA"
  ),
  group_xlim = list(
    "Harvest" = c(-3, 1.5)  
  ),
  shared_xlim_groups = c("Environmental", "Prescribed Fire", "Mastication", "Thinning"),
  xlim = c(-0.5, 0.25),  # Default for shared groups
  include_zi_effects = FALSE,
  show_component_prefix = FALSE,
  plot_height = 11,
  plot_width = 9
)

source("code/functions/marginal_effects_multimodel.R")

models <- list(
  "Breeding" = readRDS("models/pinjay_zip_breeding.rds"),
  "Non-breeding" = readRDS("models/pinjay_zip_nonbreeding.rds")
)

breeding.orig <- arrow::read_parquet("data/intermediate/pinjay_pines_ba.parquet") %>% 
  mutate(is_mobile = if_else(is_stationary, 0L, 1L)) %>% 
  filter(day_of_year >= 32 & day_of_year <= 181)
nonbreeding.orig <- arrow::read_parquet("data/intermediate/pinjay_pines_ba.parquet") %>% 
  mutate(is_mobile = if_else(is_stationary, 0L, 1L)) %>% 
  filter(day_of_year < 32 | day_of_year > 181)
original.list <- list(
  "Breeding" = breeding.orig,
  "Non-breeding" = nonbreeding.orig
)

# Compare marginal effects across models
results_main <- plot_zip_marginal_effects_multimodel(
  models = models,
  original_data_list = original.list,  # Single df used for all
  predictors = c("year_scaled", "elevation_30m_median_scaled", "pin_ba_log_scaled", 
                 "pipo_ba_log_scaled", "propwater_scaled"),
  predictor_labels = c(
    "year_scaled" = "Year", 
    "elevation_30m_median_scaled" = "Elevation (m)", 
    "pin_ba_log_scaled" = "Pinyon BA (log)", 
    "pipo_ba_log_scaled" = "Ponderosa BA (log)", 
    "propwater_scaled" = "Prop. Water"
  ),
  model_colors = c("Breeding" = "#C4956A", "Non-breeding" = "#404040"),
  n_draws = 500,
  ribbon_alpha = 0.3
)

results_main$combined_plot
ggsave("figures/main_effects_linear.png", results_main$combined_plot, width = 10, height = 12)

## let's pull out the two basal area plots
pin.pred <- filter(results_main$predictions, predictor == "pin_ba_log_scaled") |> 
  mutate(predictor_value_orig = exp(x_original) - 0.1) |> 
  select(.epred:x_original, predictor_value_orig)
pipo.pred <- filter(results_main$predictions, predictor == "pipo_ba_log_scaled") |>
  mutate(predictor_value_orig = exp(x_original) - 0.1) |> 
  select(.epred:x_original, predictor_value_orig)
library(patchwork)

## plot new
model_colors = c("Breeding" = "#C4956A", "Non-breeding" = "#404040")
p.pin <- ggplot(pin.pred, aes(x = predictor_value_orig, y = .epred,
                              color = model_id, fill = model_id)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), 
              alpha = 0.3, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = model_colors, name = "Model") +
  scale_fill_manual(values = model_colors, name = "Model") +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 0.35, base = exp(1)),
    breaks = c(0, 1, 5, 15, 50), 
    limits = c(0, 112)) +
  # scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Pinyon BA (m²/ha)",
    y = "Expected Count" 
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
    # legend.position = "none"
  )

p.pipo <- ggplot(pipo.pred, aes(x = predictor_value_orig, y = .epred,
                                color = model_id, fill = model_id)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), 
              alpha = 0.3, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = model_colors, name = "Model") +
  scale_fill_manual(values = model_colors, name = "Model") +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 0.35, base = exp(1)),
    breaks = c(0, 1, 5, 15, 50), 
    limits = c(0, 112)) +
  # scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Ponderosa BA (m²/ha)",
    y = "Expected Count" 
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

p.ba <- p.pin + p.pipo +
  plot_layout(ncol = 1, guides = "collect", axes = "collect") & theme(legend.position = "top")

ggsave("figures/basal_area_me.png", p.ba, width = 6, height = 8)

## now will add art
library(magick)

base.fig <- image_read("figures/basal_area_me.png")
pj <- image_read("figures/TerraDawson/Pimo_sketch.png") |> 
  ## resize overlay
  image_scale("16%")
pp <- image_read("figures/TerraDawson/Pipo_sketch.png") |> 
  image_scale("17%")

## combine
p <- image_composite(base.fig, pj, offset = "+1350+150") ##offset is pixel's from topleft
p <- image_composite(p, pp, offset = "+1380+1250")

## save
image_write(p, "figures/Figure4.png", quality = 100, density = 300)

## ploting marginal effects for interactions
## for setting basal area levels, let's look at distributions
hist(breeding.orig$pinyon_ba)
## add vertical lines at quantiles
abline(v = quantile(breeding.orig$pinyon_ba, probs = c(0.2, 0.5, 0.8)), col = "red")
quantile(breeding.orig$pinyon_ba, probs = c(0.2, 0.5, 0.8))

hist(breeding.orig$ponderosa_ba)
abline(v = quantile(breeding.orig$ponderosa_ba, probs = c(0.2, 0.5, 0.8)), col = "red")
quantile(breeding.orig$ponderosa_ba, probs = c(0.2, 0.5, 0.8))

hist(nonbreeding.orig$pinyon_ba)
abline(v = quantile(nonbreeding.orig$pinyon_ba, probs = c(0.2, 0.5, 0.8)), col = "red")
quantile(nonbreeding.orig$pinyon_ba, probs = c(0.2, 0.5, 0.8))

hist(nonbreeding.orig$ponderosa_ba)
abline(v = quantile(nonbreeding.orig$ponderosa_ba, probs = c(0.2, 0.5, 0.8)), col = "red")
quantile(nonbreeding.orig$ponderosa_ba, probs = c(0.2, 0.5, 0.8))

## rather than going with what is common in the dataset, let's set BA levels around what is considered under-, properly- and over-stocked (i.e., what might drive veg management decisions)

source("code/functions/marginal_effects_multimodel.R")
results.int <- plot_zip_interaction_effects_multimodel(
  models = models,
  original_data_list = original.list,
  treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
  treatment_labels = NULL,
  model_colors = c("Breeding" = "#C4956A", "Non-breeding" = "#404040"),
  yrsince_levels = c(1, 5, 10),
  pinba_levels = c(5, 35), #c(5, 15, 30),
  pipoba_levels = c(5, 35), #c(5, 15, 30),
  n_points = 51,
  ci_width = 0.9,
  n_draws = NULL,
  y_limits = c(-0.1, 8),
  x_limits = c(0, 1),
  base_font_size = 12
)

## also run with just breeding, but no harvest
int.br <- plot_zip_interaction_effects_multimodel(
  models = models[1],
  original_data_list = original.list,
  treatments = c("prescribed_fire", "mastication", "thinning"),
  treatment_labels = NULL,
  model_colors = c("Breeding" = "#C4956A", "Non-breeding" = "#404040"),
  yrsince_levels = c(1, 5, 10),
  pinba_levels = c(5, 35), #c(5, 15, 30),
  pipoba_levels = c(5, 35), #c(5, 15, 30),
  n_points = 51,
  ci_width = 0.9,
  n_draws = NULL,
  y_limits = c(-0.1, 8),
  x_limits = c(0, 1),
  base_font_size = 12
)

results.int$combined_plot_yrsince
results.int$combined_plot_pinba
results.int$combined_plot_ba

int.br$combined_plot_ba

## let's also just run for the non-breeding and harvest for that submodel
nb.har <- plot_zip_interaction_effects_multimodel(
  models = models[2],
  original_data_list = original.list[2],
  treatments = c("harvest"),
  treatment_labels = NULL,
  model_colors = c("Breeding" = "#C4956A", "Non-breeding" = "#404040"),
  yrsince_levels = c(1, 5, 10),
  pinba_levels = c(5, 35), #c(5, 15, 30),
  pipoba_levels = c(5, 35), #c(5, 15, 30),
  n_points = 51,
  ci_width = 0.9,
  n_draws = NULL,
  y_limits = c(-0.1, 8),
  x_limits = c(0, 1),
  base_font_size = 12
)

nb.har$plots_pinba

## let's tweak the combinations aesthetics a bit
library(patchwork)
p.ts <- results.int$plots_yrsince$prescribed_fire +
  results.int$plots_yrsince$mastication +
  results.int$plots_yrsince$thinning +
  plot_layout(nrow = 1, guides = "collect", axes = "collect") 

p.ba <- results.int$plots_pinba$prescribed_fire +
  results.int$plots_pinba$mastication +
  results.int$plots_pinba$thinning +
  nb.har$plots_pinba + theme(legend.position = "none") + ylab(NULL) + 
  annotate("text", x = 0.98, y = 2, angle = 90, label = "Pinyon pine") +
  results.int$plots_pipoba$prescribed_fire + labs(title = NULL) +
  results.int$plots_pipoba$mastication + labs(title = NULL) +
  results.int$plots_pipoba$thinning + labs(title = NULL) +
  results.int$plots_pipoba$harvest + labs(title = NULL) +
  annotate("text", x = 0.98, y = 2, angle = 90, label = "Ponderosa pine") +
  plot_layout(nrow = 2, guides = "collect", axes = "collect") 

ggsave("figures/seasonal_yrsince.png", p.ts, width = 10, height = 4, dpi = 300)
ggsave("figures/seasonal_ba.png", p.ba, width = 11, height = 6, dpi = 300)

## repeat for just the breeding (no harvest) version
p.ts.br <- int.br$plots_yrsince$prescribed_fire +
  int.br$plots_yrsince$mastication +
  int.br$plots_yrsince$thinning +
  plot_layout(nrow = 1, guides = "collect", axes = "collect")

p.ba.br <- int.br$plots_pinba$prescribed_fire +
  int.br$plots_pinba$mastication +
  int.br$plots_pinba$thinning +
  ## add blank space
  plot_spacer() + ylab(NULL) +
  annotate("text", x = 0.98, y = 2, angle = 90, label = "Pinyon pine") +
  int.br$plots_pipoba$prescribed_fire + labs(title = NULL) +
  int.br$plots_pipoba$mastication + labs(title = NULL) +
  int.br$plots_pipoba$thinning + labs(title = NULL) +
  plot_spacer() + ylab(NULL) +
  annotate("text", x = 0.98, y = 2, angle = 90, label = "Ponderosa pine") +
  plot_layout(nrow = 2, guides = "collect", axes = "collect")

## save
ggsave("figures/seasonal_yrsince_breeding.png", p.ts.br, width = 10, height = 4, dpi = 300)
ggsave("figures/seasonal_ba_breeding.png", p.ba.br, width = 11, height = 6, dpi = 300)

library(magick)

base.fig <- image_read("figures/seasonal_ba.png")
pj <- image_read("figures/TerraDawson/Pimo_sketch.png") |> 
  ## resize overlay
  image_scale("15%")
pp <- image_read("figures/TerraDawson/Pipo_sketch.png") |> 
  image_scale("16%")

## combine
p <- image_composite(base.fig, pj, offset = "+2800+100") ##offset is pixel's from topleft
p <- image_composite(p, pp, offset = "+2775+1200")

## save
image_write(p, "figures/Figure6.png", quality = 100, density = 300)

# ## now looking at contrasts of breeding and non-breeding
# source("code/functions/plot_zip_effect_contrasts.R")
# 
# results <- plot_zip_effect_contrasts(
#   model_a = readRDS("models/pinjay_zip_breeding.rds"),
#   model_b = readRDS("models/pinjay_zip_nonbreeding.rds"),
#   model_a_name = "Breeding",
#   model_b_name = "Non-breeding",
#   csv_output = "data/results/seasonal_contrasts_main.csv",
#   plot_output = "figures/seasonal_contrasts_main.png",
#   effects_to_include = rev(c(
#     "year_scaled",
#     "elevation_30m_median_scaled",
#     "propwater_scaled",
#     "pin_ba_log_scaled",
#     "Ipin_ba_log_scaledE2",
#     "pipo_ba_log_scaled",
#     "Ipipo_ba_log_scaledE2",
#     "prescribed_fire_bin:prescribed_fire_scaled",
#     "prescribed_fire_bin:prescribed_fire_scaled:prescribed_fire_yrsince_scaled",
#     "pin_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled",
#     "pipo_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled",
#     "mastication_bin:mastication_scaled",
#     "mastication_bin:mastication_scaled:mastication_yrsince_scaled",
#     "pin_ba_log_scaled:mastication_bin:mastication_scaled",
#     "pipo_ba_log_scaled:mastication_bin:mastication_scaled",
#     "thinning_bin:thinning_scaled",
#     "thinning_bin:thinning_scaled:thinning_yrsince_scaled",
#     "pin_ba_log_scaled:thinning_bin:thinning_scaled",
#     "pipo_ba_log_scaled:thinning_bin:thinning_scaled",
#     "harvest_bin:harvest_scaled",
#     "harvest_bin:harvest_scaled:harvest_yrsince_scaled",
#     "pin_ba_log_scaled:harvest_bin:harvest_scaled",
#     "pipo_ba_log_scaled:harvest_bin:harvest_scaled"
#   )),
#   effect_labels = c(
#     "year_scaled" = "Year",
#     "elevation_30m_median_scaled" = "Elevation",
#     "pin_ba_log_scaled" = "Pinyon Pine BA",
#     "Ipin_ba_log_scaledE2" = "Pinyon Pine BA²",
#     "pipo_ba_log_scaled" = "Ponderosa Pine BA",
#     "Ipipo_ba_log_scaledE2" = "Ponderosa Pine BA²",
#     "propwater_scaled" = "Proportion Water",
#     "prescribed_fire_bin:prescribed_fire_scaled" = "Prescribed Fire",
#     "prescribed_fire_bin:prescribed_fire_scaled:prescribed_fire_yrsince_scaled" = "Prescribed Fire x Time",
#     "pin_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled" = "Prescribed Fire x Pinyon BA",
#     "pipo_ba_log_scaled:prescribed_fire_bin:prescribed_fire_scaled" = "Prescribed Fire x Ponderosa BA",
#     "mastication_bin:mastication_scaled" = "Mastication",
#     "mastication_bin:mastication_scaled:mastication_yrsince_scaled" = "Mastication x Time",
#     "pin_ba_log_scaled:mastication_bin:mastication_scaled" = "Mastication x Pinyon BA",
#     "pipo_ba_log_scaled:mastication_bin:mastication_scaled" = "Mastication x Ponderosa BA",
#     "thinning_bin:thinning_scaled" = "Thinning",
#     "thinning_bin:thinning_scaled:thinning_yrsince_scaled" = "Thinning x Time",
#     "pin_ba_log_scaled:thinning_bin:thinning_scaled" = "Thinning x Pinyon BA",
#     "pipo_ba_log_scaled:thinning_bin:thinning_scaled" = "Thinning x Ponderosa BA",
#     "harvest_bin:harvest_scaled" = "Harvest",
#     "harvest_bin:harvest_scaled:harvest_yrsince_scaled" = "Harvest x Time",
#     "pin_ba_log_scaled:harvest_bin:harvest_scaled" = "Harvest x Pinyon BA",
#     "pipo_ba_log_scaled:harvest_bin:harvest_scaled" = "Harvest x Ponderosa BA"
#   ),
#   include_zi_effects = FALSE,
#   show_component_prefix = FALSE,
#   plot_height = 6,
#   plot_width = 12
# )
# 
# ## looking at how effects differ between breeding and non-breeding seasons
# source("code/functions/plot_seasonal_model_comparisons.R")
# 
# original_data <- arrow::read_parquet("data/intermediate/pinjay_pines_ba.parquet") %>% 
#   mutate(is_mobile = if_else(is_stationary, 0L, 1L))
# 
# results <- plot_seasonal_model_comparisons(
#   models = list(
#     "Breeding" = readRDS("models/pinjay_zip_breeding.rds"),
#     "Non-breeding" = readRDS("models/pinjay_zip_nonbreeding.rds")
#   ),
#   model_seasons = c(
#     "Breeding" = "breeding",
#     "Non-breeding" = "nonbreeding"
#   ),
#   original_data = original_data,
#   include_pipoba = TRUE,
#   treatments = c("prescribed_fire", "mastication", "thinning"),
#   prop_treated = 0.50,
#   n_draws = 500
# )
# 
# # Stacked column of years-since plots
# ggsave("figures/seasonal_yrsince.png", results$combined_plot_yrsince, 
#        width = 6, height = 10, dpi = 300)
# 
# # Two-row BA figure
# ggsave("figures/seasonal_ba_combined.png", results$combined_plot_ba,
#        width = 14, height = 8, dpi = 300)
# 
# ## looking at a version with treatment effects contrasts
# source("code/functions/plot_seasonal_treatment_effects.R")
# results <- plot_seasonal_treatment_effects(
#   models = list(
#     "Breeding" = readRDS("models/pinjay_zip_breeding.rds"),
#     "Non-breeding" = readRDS("models/pinjay_zip_nonbreeding.rds")
#   ),
#   original_data = original_data,
#   treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
#   prop_treated = 0.50,  # Compare 50% treated vs untreated
#   include_pipoba = TRUE,
#   n_draws = 500,
#   y_limits = list(yrsince = c(-2, 6), pinba = c(-2, 6), pipoba = c(-2, 6)),
# )
# 
# # Treatment effects over time
# ggsave("figures/treatment_effects_yrsince.png", results$combined_plot_yrsince,
#        width = 6, height = 10, dpi = 300)
# 
# ## currently main function fails to collect, do so here
# p <- results$plots_pinba$prescribed_fire + results$plots_pinba$mastication +
#   results$plots_pinba$thinning + results$plots_pinba$harvest +
#   results$plots_pipoba$prescribed_fire + results$plots_pipoba$mastication +
#   results$plots_pipoba$thinning + results$plots_pipoba$harvest +
#   plot_layout(guides = "collect", axes = "collect", nrow = 2)
# 
# # Two-row BA figure
# ggsave("figures/treatment_effects_ba.png", p,
#        width = 14, height = 8, dpi = 300)

#### Spotted Owl ####

## prep spotted owl data
source("code/functions/sp_data_org_generalized.R")
sp_data_org("spoowl")


#### Black-backed woodpecker ####

## prep black-backed woodpecker data
source("code/functions/sp_data_org_generalized.R")
sp_data_org("bkbwoo")


#### Clark's Nutcracker ####
## prep clark's nutcracker data
source("code/functions/sp_data_org_generalized.R")
sp_data_org("clanut")
