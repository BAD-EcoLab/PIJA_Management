## Purpose: Pinyon Jay data management and analysis dashboard
## Project: PIJA_Management

## plot the prevalence of different treatment types with respect to pinyon and ponderosa pine density
source("code/functions/pine_treatments_plot.R")
pine_treatments_plot()


## Run models for the breeding and non-breeding seasons
## These take bout 20 hrs to run (10 each)
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
    "Breeding" = readRDS("models/pinjay_zip_breeding.rds"),
    "Non-breeding" = readRDS("models/pinjay_zip_nonbreeding.rds")
  ),
  csv_output = "data/seasonal_comparison_main.csv",
  plot_output = "figures/Figure_3.png",
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

breeding.orig <- arrow::read_parquet("data/model_data.parquet") %>% 
  mutate(is_mobile = if_else(is_stationary, 0L, 1L)) %>% 
  filter(day_of_year >= 32 & day_of_year <= 181)
nonbreeding.orig <- arrow::read_parquet("data/model_data.parquet") %>% 
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
  labs(
    x = "Pinyon BA (m²/ha)",
    y = "Expected Count" 
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
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

## note that this lacks the artwork used in the publication
ggsave("figures/Figure_4.png", p.ba, width = 6, height = 8)


## ploting marginal effects for interactions
## setting BA levels around what is considered under-, and over-stocked (i.e., what might drive veg management decisions)

source("code/functions/marginal_effects_multimodel.R")
results.int <- plot_zip_interaction_effects_multimodel(
  models = models,
  original_data_list = original.list,
  treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
  treatment_labels = NULL,
  model_colors = c("Breeding" = "#C4956A", "Non-breeding" = "#404040"),
  yrsince_levels = c(1, 5, 10),
  pinba_levels = c(5, 35), 
  pipoba_levels = c(5, 35), 
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
  pinba_levels = c(5, 35), 
  pipoba_levels = c(5, 35), 
  n_points = 51,
  ci_width = 0.9,
  n_draws = NULL,
  y_limits = c(-0.1, 8),
  x_limits = c(0, 1),
  base_font_size = 12
)

## let's also just run for the non-breeding and harvest for that submodel
nb.har <- plot_zip_interaction_effects_multimodel(
  models = models[2],
  original_data_list = original.list[2],
  treatments = c("harvest"),
  treatment_labels = NULL,
  model_colors = c("Breeding" = "#C4956A", "Non-breeding" = "#404040"),
  yrsince_levels = c(1, 5, 10),
  pinba_levels = c(5, 35), 
  pipoba_levels = c(5, 35), 
  n_points = 51,
  ci_width = 0.9,
  n_draws = NULL,
  y_limits = c(-0.1, 8),
  x_limits = c(0, 1),
  base_font_size = 12
)

## Adjusting combinations aesthetics a bit
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

ggsave("figures/Figure_5.png", p.ts, width = 10, height = 4, dpi = 300)
ggsave("figures/Figure_6.png", p.ba, width = 11, height = 6, dpi = 300)


