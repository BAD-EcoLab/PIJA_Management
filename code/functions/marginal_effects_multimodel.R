## Multi-model marginal effects plots for ZIP models

library(tidyverse)
library(brms)
library(tidybayes)
library(ggdist)
library(modelr)
library(patchwork)
library(scales)

# Extract scaling parameters for back-transformation
get_scaling_params <- function(original_data, #Original UNSCALED data
                               vars_to_scale #Character vector of variable names (without _scaled suffix)
                               ) {
  params <- list()
  for (var in vars_to_scale) {
    if (var %in% names(original_data)) {
      params[[paste0(var, "_scaled")]] <- list(
        mean = mean(original_data[[var]], na.rm = TRUE),
        sd = sd(original_data[[var]], na.rm = TRUE),
        orig_name = var
      )
    }
  }
  return(params)
}



# Plot marginal effects for non-interacting predictors across multiple models
plot_zip_marginal_effects_multimodel <- function(
    models,
    original_data_list,
    predictors = c("year_scaled", "day_of_year_scaled", "effort_hrs_scaled",
                   "num_observers_scaled", "cci_scaled", 
                   "elevation_30m_median_scaled", "propwater_scaled",
                   "pin_ba_log_scaled"),
    predictor_labels = NULL,
    model_colors = NULL,
    model_linetypes = NULL,
    n_points = 51,
    ci_width = 0.9,
    n_draws = NULL,
    ncol = 2,
    show_ribbon = TRUE,
    ribbon_alpha = 0.15,
    base_font_size = 11
) {
  
  # Validate inputs
  if (!is.list(models) || is.null(names(models))) {
    stop("'models' must be a named list, e.g., list('Full Year' = m1, 'Breeding' = m2)")
  }
  model_names <- names(models)
  message("Processing ", length(models), " models: ", paste(model_names, collapse = ", "))
  
  # Handle original_data_list - allow single df or named list
  if (is.data.frame(original_data_list)) {
    # Single df provided - use for all models
    original_data_list <- setNames(
      replicate(length(models), original_data_list, simplify = FALSE),
      model_names
    )
  } else if (!is.list(original_data_list) || is.null(names(original_data_list))) {
    stop("'original_data_list' must be a named list matching 'models', or a single data frame")
  }
  
  # Set up colors
  if (is.null(model_colors)) {
    default_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
    model_colors <- setNames(default_colors[1:length(model_names)], model_names)
  }
  
  # Set up linetypes (optional)
  if (is.null(model_linetypes)) {
    model_linetypes <- setNames(rep("solid", length(model_names)), model_names)
  }
  
  # Default labels
  if (is.null(predictor_labels)) {
    predictor_labels <- setNames(
      c("Year", "Day of Year", "Effort (hrs)", "# Observers", 
        "CCI", "Elevation (m)", "Prop. Water", 
        "Pinyon BA (log)", "Ponderosa Pine BA (log)"),
      c("year_scaled", "day_of_year_scaled", "effort_hrs_scaled",
        "num_observers_scaled", "cci_scaled", 
        "elevation_30m_median_scaled", "propwater_scaled",
        "pin_ba_log_scaled", "pipo_ba_log_scaled")
    )
  }
  
  # Binary treatment indicators - set to 0 for baseline
  binary_vars <- c("prescribed_fire_bin", "mastication_bin", "thinning_bin",
                   "harvest_bin", "clearcut_bin", "is_mobile")
  
  # Storage for all predictions
  all_preds <- list()
  
  # Process each model
  for (model_name in model_names) {
    message("  Processing model: ", model_name)
    
    model <- models[[model_name]]
    original_data <- original_data_list[[model_name]]
    model_data <- model$data
    
    # Get all scaled vars and create baseline
    all_scaled_vars <- names(model_data)[grepl("_scaled$", names(model_data))]
    
    baseline <- model_data %>%
      summarise(across(all_of(all_scaled_vars), mean, na.rm = TRUE))
    
    for (bv in binary_vars) {
      if (bv %in% names(model_data)) {
        baseline[[bv]] <- 0
      }
    }
    
    # Get scaling params for this model's original data
    vars_for_scaling <- gsub("_scaled$", "", predictors)
    scaling_params <- get_scaling_params(original_data, vars_for_scaling)
    
    # Process each predictor
    for (pred in predictors) {
      if (!pred %in% names(model_data)) {
        message("    Skipping ", pred, " - not in model")
        next
      }
      
      pred_range <- range(model_data[[pred]], na.rm = TRUE)
      
      pred_grid <- baseline %>%
        slice(rep(1, n_points)) %>%
        mutate(!!pred := seq(pred_range[1], pred_range[2], length.out = n_points))
      
      pred_grid$ecoregion <- "new"
      
      # Extract predictions
      preds <- pred_grid %>%
        add_epred_draws(
          model,
          re_formula = NA,
          allow_new_levels = TRUE,
          ndraws = n_draws
        )
      
      # Summarize
      pred_summary <- preds %>%
        group_by(across(all_of(pred))) %>%
        median_qi(.epred, 
                  .width = ci_width) %>% 
        ungroup() %>%
        mutate(
          predictor = pred,
          model_id = model_name
        )
      
      # Back-transform to original scale
      if (pred %in% names(scaling_params)) {
        sp <- scaling_params[[pred]]
        pred_summary <- pred_summary %>%
          mutate(x_original = .data[[pred]] * sp$sd + sp$mean)
      } else {
        pred_summary <- pred_summary %>%
          mutate(x_original = .data[[pred]])
      }
      
      all_preds[[paste0(model_name, "_", pred)]] <- pred_summary
    }
  }
  
  # Combine all predictions
  all_preds_df <- bind_rows(all_preds) %>%
    mutate(model_id = factor(model_id, levels = model_names))
  
  message("Creating plots...")
  
  # Create individual plots
  plots <- list()
  
  for (pred in predictors) {
    pred_data <- all_preds_df %>%
      filter(predictor == pred)
    
    if (nrow(pred_data) == 0) next
    
    x_label <- if (pred %in% names(predictor_labels)) {
      predictor_labels[[pred]]
    } else {
      gsub("_scaled$", "", pred)
    }
    
    pred_index <- which(predictors == pred)
    is_left_column <- ((pred_index - 1) %% ncol) == 0
    
    p <- ggplot(pred_data, aes(x = x_original, y = .epred, 
                               color = model_id, fill = model_id,
                               linetype = model_id))
    
    if (show_ribbon) {
      p <- p + geom_ribbon(aes(ymin = .lower, ymax = .upper), 
                           alpha = ribbon_alpha, color = NA)
    }
    
    p <- p +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors, name = "Model") +
      scale_fill_manual(values = model_colors, name = "Model") +
      scale_linetype_manual(values = model_linetypes, name = "Model") +
      scale_x_continuous(labels = label_number(accuracy = 1)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      labs(
        x = x_label,
        y = if (is_left_column) "Expected Count" else NULL
      ) +
      theme_bw(base_size = base_font_size) +
      theme(
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.5, "cm")
        # legend.position = "none"
      )
    
    plots[[pred]] <- p
  }
  
  # Create combined plot
  combined <- wrap_plots(plots, ncol = ncol, guides = "collect") +
    plot_annotation(
      title = "Marginal Effects: Non-interacting Predictors",
      subtitle = paste0("Median with ", ci_width * 100, "% CI | ", 
                        paste(model_names, collapse = " vs "))
    ) &
    theme(legend.position = "bottom")
  
  return(list(
    predictions = all_preds_df,
    plots = plots,
    combined_plot = combined
  ))
}


# Plot treatment interaction effects across multiple models (non-flipped)

# X-axis: Proportion treated
# Color/lines: Different years-since or BA levels
# Panels: Different models shown together via color/linetype
# Includes credible intervals via geom_ribbon
plot_zip_interaction_effects_multimodel <- function(
    models,
    original_data_list,
    treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
    treatment_labels = NULL,
    model_colors = NULL,
    yrsince_levels = c(1, 5, 10),
    pinba_levels = c(1, 5, 10),
    pipoba_levels = c(1, 5, 10),
    n_points = 51,
    ci_width = 0.9,
    n_draws = NULL,
    y_limits = c(-0.1, 10),
    x_limits = c(0, 1),
    base_font_size = 10
) {
  
  # Validate inputs
  if (!is.list(models) || is.null(names(models))) {
    stop("'models' must be a named list")
  }
  model_names <- names(models)
  message("Processing ", length(models), " models: ", paste(model_names, collapse = ", "))
  
  # Handle original_data_list
  if (is.data.frame(original_data_list)) {
    original_data_list <- setNames(
      replicate(length(models), original_data_list, simplify = FALSE),
      model_names
    )
  }
  
  # Set up colors
  if (is.null(model_colors)) {
    default_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
    model_colors <- setNames(default_colors[1:length(model_names)], model_names)
  }
  
  # Default treatment labels
  if (is.null(treatment_labels)) {
    treatment_labels <- c(
      "prescribed_fire" = "Prescribed Fire",
      "mastication" = "Mastication", 
      "thinning" = "Thinning",
      "harvest" = "Harvest",
      "clearcut" = "Clearcut"
    )
  }
  
  # Validate basal area levels
  if (any(pinba_levels <= -0.1)) {
    stop("pinba_levels must be greater than -0.1 to avoid log of negative numbers")
  }
  if (any(pipoba_levels <= -0.1)) {
    stop("pipoba_levels must be greater than -0.1 to avoid log of negative numbers")
  }
  if (length(pinba_levels) > 6) {
    warning("pinba_levels has more than 6 values; only first 6 will be distinguishable by linetype")
    pinba_levels <- pinba_levels[1:6]
  }
  if (length(pipoba_levels) > 6) {
    warning("pipoba_levels has more than 6 values; only first 6 will be distinguishable by linetype")
    pipoba_levels <- pipoba_levels[1:6]
  }
  
  binary_vars <- c("prescribed_fire_bin", "mastication_bin", "thinning_bin",
                   "harvest_bin", "clearcut_bin", "is_mobile")
  
  # Convert absolute BA values to log scale (with 0.1 offset to avoid log(0))
  pinba_values_log <- log(pinba_levels + 0.1)
  pipoba_values_log <- log(pipoba_levels + 0.1)
  
  # Create labels using absolute BA values
  pinba_labels <- as.character(pinba_levels)
  pipoba_labels <- as.character(pipoba_levels)
  
  # Storage
  all_preds <- list()
  plots_yrsince <- list()
  plots_pinba <- list()
  plots_pipoba <- list()
  
  for (trt in treatments) {
    message("Processing treatment: ", trt)
    
    trt_scaled <- paste0(trt, "_scaled")
    trt_bin <- paste0(trt, "_bin")
    yrsince_scaled <- paste0(trt, "_yrsince_scaled")
    
    trt_label <- ifelse(trt %in% names(treatment_labels),
                        treatment_labels[trt], trt)
    
    # Collect predictions across models for this treatment
    preds_yrsince_all <- list()
    preds_pinba_all <- list()
    preds_pipoba_all <- list()
    
    for (model_name in model_names) {
      model <- models[[model_name]]
      original_data <- original_data_list[[model_name]]
      model_data <- model$data
      
      # Skip if treatment not in model
      if (!trt_scaled %in% names(model_data)) {
        message("  Skipping ", trt, " for model ", model_name, " - not found")
        next
      }
      
      # Get scaling params
      all_scaled_vars <- names(model_data)[grepl("_scaled$", names(model_data))]
      
      baseline <- model_data %>%
        summarise(across(all_of(all_scaled_vars), mean, na.rm = TRUE))
      
      for (bv in binary_vars) {
        if (bv %in% names(model_data)) baseline[[bv]] <- 0
      }
      
      # Scaling parameters for this model
      sp_trt <- list(
        mean = mean(original_data[[trt]], na.rm = TRUE),
        sd = sd(original_data[[trt]], na.rm = TRUE)
      )
      
      yrsince_var <- paste0(trt, "_yrsince")
      sp_yrsince <- list(
        mean = mean(original_data[[yrsince_var]], na.rm = TRUE),
        sd = sd(original_data[[yrsince_var]], na.rm = TRUE)
      )
      
      sp_pinba <- list(
        mean = mean(original_data$pin_ba_log, na.rm = TRUE),
        sd = sd(original_data$pin_ba_log, na.rm = TRUE)
      )
      # Use log-transformed absolute BA values instead of quantiles
      pinba_values_orig <- pinba_values_log
      
      sp_pipoba <- list(
        mean = mean(original_data$pipo_ba_log, na.rm = TRUE),
        sd = sd(original_data$pipo_ba_log, na.rm = TRUE)
      )
      # Use log-transformed absolute BA values instead of quantiles
      pipoba_values_orig <- pipoba_values_log
      
      trt_range <- range(model_data[[trt_scaled]], na.rm = TRUE)
      
      # Years-since interaction
      yrsince_scaled_values <- (yrsince_levels - sp_yrsince$mean) / sp_yrsince$sd
      
      preds_yrs <- map_dfr(seq_along(yrsince_levels), function(i) {
        pred_grid <- baseline %>%
          slice(rep(1, n_points)) %>%
          mutate(
            !!trt_scaled := seq(trt_range[1], trt_range[2], length.out = n_points),
            !!trt_bin := 1,
            !!yrsince_scaled := yrsince_scaled_values[i]
          )
        pred_grid$ecoregion <- "new"
        
        preds <- pred_grid %>%
          add_epred_draws(model, re_formula = NA, allow_new_levels = TRUE, ndraws = n_draws) %>%
          group_by(across(all_of(trt_scaled))) %>%
          median_qi(.epred, .width = ci_width) %>%
          ungroup() %>%
          mutate(
            yrsince_raw = yrsince_levels[i],
            yrsince_label = paste0(yrsince_levels[i], " yr")
          )
        return(preds)
      })
      
      preds_yrs <- preds_yrs %>%
        mutate(
          x_original = .data[[trt_scaled]] * sp_trt$sd + sp_trt$mean,
          treatment = trt,
          model_id = model_name,
          yrsince_label = factor(yrsince_label, levels = paste0(sort(yrsince_levels), " yr"))
        )
      
      preds_yrsince_all[[model_name]] <- preds_yrs
      
      # Pinyon BA interaction
      pinba_scaled_values <- (pinba_values_orig - sp_pinba$mean) / sp_pinba$sd
      
      preds_ba <- map_dfr(seq_along(pinba_levels), function(i) {
        pred_grid <- baseline %>%
          slice(rep(1, n_points)) %>%
          mutate(
            !!trt_scaled := seq(trt_range[1], trt_range[2], length.out = n_points),
            !!trt_bin := 1,
            pin_ba_log_scaled = pinba_scaled_values[i]
          )
        pred_grid$ecoregion <- "new"
        
        preds <- pred_grid %>%
          add_epred_draws(model, re_formula = NA, allow_new_levels = TRUE, ndraws = n_draws) %>%
          group_by(across(all_of(trt_scaled))) %>%
          median_qi(.epred, .width = ci_width) %>%
          ungroup() %>%
          mutate(pinba_label = pinba_labels[i])
        return(preds)
      })
      
      preds_ba <- preds_ba %>%
        mutate(
          x_original = .data[[trt_scaled]] * sp_trt$sd + sp_trt$mean,
          treatment = trt,
          model_id = model_name,
          pinba_label = factor(pinba_label, levels = pinba_labels)
        )
      
      preds_pinba_all[[model_name]] <- preds_ba
      
      # Ponderosa Pine BA interaction
      pipoba_scaled_values <- (pipoba_values_orig - sp_pipoba$mean) / sp_pipoba$sd
      
      preds_pipoba <- map_dfr(seq_along(pipoba_levels), function(i) {
        pred_grid <- baseline %>%
          slice(rep(1, n_points)) %>%
          mutate(
            !!trt_scaled := seq(trt_range[1], trt_range[2], length.out = n_points),
            !!trt_bin := 1,
            pipo_ba_log_scaled = pipoba_scaled_values[i]
          )
        pred_grid$ecoregion <- "new"
        
        preds <- pred_grid %>%
          add_epred_draws(model, re_formula = NA, allow_new_levels = TRUE, ndraws = n_draws) %>%
          group_by(across(all_of(trt_scaled))) %>%
          median_qi(.epred, .width = ci_width) %>%
          ungroup() %>%
          mutate(pipoba_label = pipoba_labels[i])
        return(preds)
      })
      
      preds_pipoba <- preds_pipoba %>%
        mutate(
          x_original = .data[[trt_scaled]] * sp_trt$sd + sp_trt$mean,
          treatment = trt,
          model_id = model_name,
          pipoba_label = factor(pipoba_label, levels = pipoba_labels)
        )
      
      preds_pipoba_all[[model_name]] <- preds_pipoba
    }
    
    # Combine across models
    preds_yrsince_combined <- bind_rows(preds_yrsince_all) %>%
      mutate(model_id = factor(model_id, levels = model_names))
    
    preds_pinba_combined <- bind_rows(preds_pinba_all) %>%
      mutate(model_id = factor(model_id, levels = model_names))
    
    preds_pipoba_combined <- bind_rows(preds_pipoba_all) %>%
      mutate(model_id = factor(model_id, levels = model_names))
    
    all_preds[[paste0(trt, "_yrsince")]] <- preds_yrsince_combined
    all_preds[[paste0(trt, "_pinba")]] <- preds_pinba_combined
    all_preds[[paste0(trt, "_pipoba")]] <- preds_pipoba_combined
    
    # Determine position
    trt_index <- which(treatments == trt)
    is_first_col <- trt_index == 1
    is_last_col <- trt_index == length(treatments)
    
    # Create years-since plot
    p_yrs <- ggplot(preds_yrsince_combined, 
                    aes(x = x_original, y = .epred, 
                        color = model_id, linetype = yrsince_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = model_id), 
                  alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors, name = "Model") +
      scale_fill_manual(values = model_colors, guide = "none") +
      scale_linetype_manual(
        values = c("1 yr" = "solid", "5 yr" = "dashed", "10 yr" = "dotted"),
        name = "Years Since"
      ) +
      scale_y_continuous(
        trans = pseudo_log_trans(sigma = 0.5),
        breaks = c(0, 1, 3, 6, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = x_limits, ylim = y_limits) +
      scale_x_continuous(labels = function(x) x * 100,  # Multiply values by 100
                         breaks = seq(0, 1, 0.25)) +      # Define break points
      labs(
        x = "Landscape Percent Treated",
        y = if (is_first_col) "Expected Count" else NULL,
        title = trt_label
      ) +
      theme_bw(base_size = base_font_size) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        legend.key.width = unit(1.5, "cm")
      )
    
    plots_yrsince[[trt]] <- p_yrs
    
    # Create pinyon BA plot
    p_pinba <- ggplot(preds_pinba_combined,
                      aes(x = x_original, y = .epred,
                          color = model_id, linetype = pinba_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = model_id),
                  alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors, name = "Model") +
      scale_fill_manual(values = model_colors, guide = "none") +
      scale_linetype_manual(
        values = setNames(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")[seq_along(pinba_labels)], pinba_labels),
        name = "Basal Area\n(m²/ha)"
      ) +
      scale_y_continuous(
        trans = pseudo_log_trans(sigma = 0.5),
        breaks = c(0, 1, 3, 6, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = x_limits, ylim = y_limits) +
      scale_x_continuous(labels = function(x) x * 100,  # Multiply values by 100
                         breaks = seq(0, 1, 0.25)) +   # Define break points
      labs(
        x = "Landscape Percent Treated",
        y = if (is_first_col) "Expected Count" else NULL,
        title = trt_label
      ) +
      theme_bw(base_size = base_font_size) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        legend.key.width = unit(1.5, "cm")
      )
    
    plots_pinba[[trt]] <- p_pinba
    
    # Create ponderosa pine BA plot
    p_pipoba <- ggplot(preds_pipoba_combined,
                       aes(x = x_original, y = .epred,
                           color = model_id, linetype = pipoba_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = model_id),
                  alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors, name = "Model") +
      scale_fill_manual(values = model_colors, guide = "none") +
      scale_linetype_manual(
        values = setNames(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")[seq_along(pipoba_labels)], pipoba_labels),
        name = "Basal Area\n(m²/ha)"
      ) +
      scale_y_continuous(
        trans = pseudo_log_trans(sigma = 0.5),
        breaks = c(0, 1, 3, 6, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = x_limits, ylim = y_limits) +
      scale_x_continuous(labels = function(x) x * 100,  # Multiply values by 100
                         breaks = seq(0, 1, 0.25)) +     # Define break points
      labs(
        x = "Landscape Percent Treated",
        y = if (is_first_col) "Expected Count" else NULL,
        title = trt_label
      ) +
      theme_bw(base_size = base_font_size) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        legend.key.width = unit(1.5, "cm")
      )
    
    plots_pipoba[[trt]] <- p_pipoba
  }
  
  # Combine predictions
  all_preds_df <- bind_rows(all_preds)
  
  # Combined plots
  combined_yrsince <- wrap_plots(plots_yrsince, nrow = 1, guides = "collect")
  combined_pinba <- wrap_plots(plots_pinba, nrow = 1, guides = "collect")
  combined_pipoba <- wrap_plots(plots_pipoba, nrow = 1, guides = "collect")
  
  # Combined BA plot with both pine species (2 rows)
  combined_ba <- wrap_plots(
    combined_pinba,
    combined_pipoba,
    nrow = 2,
    guides = "collect",
    axes = "collect"
  ) +
    plot_annotation(
      title = "Basal Area Interactions",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
    )
  
  return(list(
    predictions = all_preds_df,
    plots_yrsince = plots_yrsince,
    plots_pinba = plots_pinba,
    plots_pipoba = plots_pipoba,
    combined_plot_yrsince = combined_yrsince,
    combined_plot_pinba = combined_pinba,
    combined_plot_pipoba = combined_pipoba,
    combined_plot_ba = combined_ba
  ))
}


# Plot treatment interaction effects across multiple models (flipped)

# X-axis: Years since treatment or Basal Area
# Color/lines: Different proportion-treated levels
# Faceted or grouped by model
# Includes credible intervals via geom_ribbon
plot_zip_interaction_effects_flipped_multimodel <- function(
    models,
    original_data_list,
    treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
    treatment_labels = NULL,
    model_colors = NULL,
    prop_treated_levels = c(0.25, 0.50, 0.75),
    yrsince_range = c(1, 10),
    n_points = 51,
    ci_width = 0.9,
    n_draws = NULL,
    filter_to_observed = TRUE,
    y_limits = c(-0.1, 6),
    ba_x_limits = c(0, 100),
    base_font_size = 10
) {
  
  # Validate inputs
  if (!is.list(models) || is.null(names(models))) {
    stop("'models' must be a named list")
  }
  model_names <- names(models)
  message("Processing ", length(models), " models: ", paste(model_names, collapse = ", "))
  
  # Handle original_data_list
  if (is.data.frame(original_data_list)) {
    original_data_list <- setNames(
      replicate(length(models), original_data_list, simplify = FALSE),
      model_names
    )
  }
  
  # Set up colors
  if (is.null(model_colors)) {
    default_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
    model_colors <- setNames(default_colors[1:length(model_names)], model_names)
  }
  
  # Default treatment labels
  if (is.null(treatment_labels)) {
    treatment_labels <- c(
      "prescribed_fire" = "Prescribed Fire",
      "mastication" = "Mastication", 
      "thinning" = "Thinning",
      "harvest" = "Harvest",
      "clearcut" = "Clearcut"
    )
  }
  
  # Proportion labels
  prop_labels <- paste0(prop_treated_levels * 100, "%")
  
  binary_vars <- c("prescribed_fire_bin", "mastication_bin", "thinning_bin",
                   "harvest_bin", "clearcut_bin", "is_mobile")
  
  # BA breaks for pseudo-log scale
  ba_breaks_orig <- c(0, 1, 5, 25, 100)
  ba_breaks_labels <- c("0", "1", "5", "25", "100")
  
  # Storage
  all_preds <- list()
  plots_yrsince <- list()
  plots_pinba <- list()
  plots_pipoba <- list()
  
  for (trt in treatments) {
    message("Processing treatment: ", trt)
    
    trt_scaled <- paste0(trt, "_scaled")
    trt_bin <- paste0(trt, "_bin")
    yrsince_scaled <- paste0(trt, "_yrsince_scaled")
    
    trt_label <- ifelse(trt %in% names(treatment_labels),
                        treatment_labels[trt], trt)
    
    # Collect predictions across models
    preds_yrsince_all <- list()
    preds_pinba_all <- list()
    preds_pipoba_all <- list()
    
    for (model_name in model_names) {
      model <- models[[model_name]]
      original_data <- original_data_list[[model_name]]
      model_data <- model$data
      
      # Skip if treatment not in model
      if (!trt_scaled %in% names(model_data)) {
        message("  Skipping ", trt, " for model ", model_name, " - not found")
        next
      }
      
      # Get baseline
      all_scaled_vars <- names(model_data)[grepl("_scaled$", names(model_data))]
      
      baseline <- model_data %>%
        summarise(across(all_of(all_scaled_vars), mean, na.rm = TRUE))
      
      for (bv in binary_vars) {
        if (bv %in% names(model_data)) baseline[[bv]] <- 0
      }
      
      # Scaling parameters
      sp_trt <- list(
        mean = mean(original_data[[trt]], na.rm = TRUE),
        sd = sd(original_data[[trt]], na.rm = TRUE)
      )
      
      yrsince_var <- paste0(trt, "_yrsince")
      sp_yrsince <- list(
        mean = mean(original_data[[yrsince_var]], na.rm = TRUE),
        sd = sd(original_data[[yrsince_var]], na.rm = TRUE)
      )
      
      sp_pinba <- list(
        mean = mean(original_data$pin_ba_log, na.rm = TRUE),
        sd = sd(original_data$pin_ba_log, na.rm = TRUE)
      )
      
      sp_pipoba <- list(
        mean = mean(original_data$pipo_ba_log, na.rm = TRUE),
        sd = sd(original_data$pipo_ba_log, na.rm = TRUE)
      )
      
      # Filter prop levels to observed range
      max_prop_observed <- max(original_data[[trt]], na.rm = TRUE)
      
      if (filter_to_observed) {
        valid_prop_idx <- which(prop_treated_levels <= max_prop_observed)
        if (length(valid_prop_idx) == 0) {
          message("  Skipping ", trt, " for ", model_name, ": no prop levels in observed range (max = ",
                  round(max_prop_observed, 2), ")")
          next
        }
        trt_prop_levels <- prop_treated_levels[valid_prop_idx]
        trt_prop_labels <- paste0(trt_prop_levels * 100, "%")
        message("    Using prop levels: ", paste(trt_prop_labels, collapse = ", "))
      } else {
        trt_prop_levels <- prop_treated_levels
        trt_prop_labels <- prop_labels
      }
      
      trt_prop_scaled <- (trt_prop_levels - sp_trt$mean) / sp_trt$sd
      
      # Pinyon BA range
      pinba_range_orig <- c(min(original_data$pin_ba_log, na.rm = TRUE),
                            max(original_data$pin_ba_log, na.rm = TRUE))
      
      # Ponderosa BA range
      pipoba_range_orig <- c(min(original_data$pipo_ba_log, na.rm = TRUE),
                             max(original_data$pipo_ba_log, na.rm = TRUE))
      
      # Years-since on x-axis
      yrsince_range_scaled <- (yrsince_range - sp_yrsince$mean) / sp_yrsince$sd
      
      preds_yrs <- map_dfr(seq_along(trt_prop_levels), function(i) {
        pred_grid <- baseline %>%
          slice(rep(1, n_points)) %>%
          mutate(
            !!yrsince_scaled := seq(yrsince_range_scaled[1], yrsince_range_scaled[2], length.out = n_points),
            !!trt_scaled := trt_prop_scaled[i],
            !!trt_bin := 1
          )
        pred_grid$ecoregion <- "new"
        
        preds <- pred_grid %>%
          add_epred_draws(model, re_formula = NA, allow_new_levels = TRUE, ndraws = n_draws) %>%
          group_by(across(all_of(yrsince_scaled))) %>%
          median_qi(.epred, .width = ci_width) %>%
          ungroup() %>%
          mutate(
            prop_raw = trt_prop_levels[i],
            prop_label = trt_prop_labels[i]
          )
        return(preds)
      })
      
      preds_yrs <- preds_yrs %>%
        mutate(
          x_original = .data[[yrsince_scaled]] * sp_yrsince$sd + sp_yrsince$mean,
          treatment = trt,
          model_id = model_name,
          prop_label = factor(prop_label, levels = prop_labels)
        )
      
      preds_yrsince_all[[model_name]] <- preds_yrs
      
      # Pinyon BA on x-axis
      pinba_range_scaled <- (pinba_range_orig - sp_pinba$mean) / sp_pinba$sd
      
      preds_ba <- map_dfr(seq_along(trt_prop_levels), function(i) {
        pred_grid <- baseline %>%
          slice(rep(1, n_points)) %>%
          mutate(
            pin_ba_log_scaled = seq(pinba_range_scaled[1], pinba_range_scaled[2], length.out = n_points),
            !!trt_scaled := trt_prop_scaled[i],
            !!trt_bin := 1
          )
        pred_grid$ecoregion <- "new"
        
        preds <- pred_grid %>%
          add_epred_draws(model, re_formula = NA, allow_new_levels = TRUE, ndraws = n_draws) %>%
          group_by(pin_ba_log_scaled) %>%
          median_qi(.epred, .width = ci_width) %>%
          ungroup() %>%
          mutate(
            prop_raw = trt_prop_levels[i],
            prop_label = trt_prop_labels[i]
          )
        return(preds)
      })
      
      preds_ba <- preds_ba %>%
        mutate(
          x_original = exp(pin_ba_log_scaled * sp_pinba$sd + sp_pinba$mean),
          treatment = trt,
          model_id = model_name,
          prop_label = factor(prop_label, levels = prop_labels)
        )
      
      preds_pinba_all[[model_name]] <- preds_ba
      
      # Ponderosa Pine BA on x-axis
      pipoba_range_scaled <- (pipoba_range_orig - sp_pipoba$mean) / sp_pipoba$sd
      
      preds_pipoba <- map_dfr(seq_along(trt_prop_levels), function(i) {
        pred_grid <- baseline %>%
          slice(rep(1, n_points)) %>%
          mutate(
            pipo_ba_log_scaled = seq(pipoba_range_scaled[1], pipoba_range_scaled[2], length.out = n_points),
            !!trt_scaled := trt_prop_scaled[i],
            !!trt_bin := 1
          )
        pred_grid$ecoregion <- "new"
        
        preds <- pred_grid %>%
          add_epred_draws(model, re_formula = NA, allow_new_levels = TRUE, ndraws = n_draws) %>%
          group_by(pipo_ba_log_scaled) %>%
          median_qi(.epred, .width = ci_width) %>%
          ungroup() %>%
          mutate(
            prop_raw = trt_prop_levels[i],
            prop_label = trt_prop_labels[i]
          )
        return(preds)
      })
      
      preds_pipoba <- preds_pipoba %>%
        mutate(
          x_original = exp(pipo_ba_log_scaled * sp_pipoba$sd + sp_pipoba$mean),
          treatment = trt,
          model_id = model_name,
          prop_label = factor(prop_label, levels = prop_labels)
        )
      
      preds_pipoba_all[[model_name]] <- preds_pipoba
    }
    
    # Combine across models
    preds_yrsince_combined <- bind_rows(preds_yrsince_all) %>%
      mutate(model_id = factor(model_id, levels = model_names))
    
    preds_pinba_combined <- bind_rows(preds_pinba_all) %>%
      mutate(model_id = factor(model_id, levels = model_names))
    
    preds_pipoba_combined <- bind_rows(preds_pipoba_all) %>%
      mutate(model_id = factor(model_id, levels = model_names))
    
    all_preds[[paste0(trt, "_yrsince_flipped")]] <- preds_yrsince_combined
    all_preds[[paste0(trt, "_pinba_flipped")]] <- preds_pinba_combined
    all_preds[[paste0(trt, "_pipoba_flipped")]] <- preds_pipoba_combined
    
    # Determine position
    trt_index <- which(treatments == trt)
    is_first_col <- trt_index == 1
    is_last_col <- trt_index == length(treatments)
    
    # Create years-since plot (flipped)
    p_yrs <- ggplot(preds_yrsince_combined,
                    aes(x = x_original, y = .epred,
                        color = model_id, linetype = prop_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = model_id),
                  alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors, name = "Model") +
      scale_fill_manual(values = model_colors, guide = "none") +
      scale_linetype_manual(
        values = setNames(
          c("solid", "dashed", "dotted", "longdash")[1:length(prop_labels)],
          prop_labels
        ),
        name = "% Treated",
        drop = FALSE
      ) +
      scale_x_continuous(labels = label_number(accuracy = 1)) +
      scale_y_continuous(
        trans = pseudo_log_trans(sigma = 1),
        breaks = c(0, 1, 2, 4, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(ylim = y_limits) +
      labs(
        x = "Years Since Treatment",
        y = if (is_first_col) "Expected Count" else NULL,
        title = trt_label
      ) +
      theme_bw(base_size = base_font_size) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = if (is_last_col) "right" else "none",
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        legend.key.width = unit(1.5, "cm")
      )
    
    plots_yrsince[[trt]] <- p_yrs
    
    # Create pinyon BA plot (flipped)
    p_pinba <- ggplot(preds_pinba_combined,
                      aes(x = x_original, y = .epred,
                          color = model_id, linetype = prop_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = model_id),
                  alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors, name = "Model") +
      scale_fill_manual(values = model_colors, guide = "none") +
      scale_linetype_manual(
        values = setNames(
          c("solid", "dashed", "dotted", "longdash")[1:length(prop_labels)],
          prop_labels
        ),
        name = "% Treated",
        drop = FALSE
      ) +
      scale_x_continuous(
        trans = pseudo_log_trans(sigma = 0.5, base = exp(1)),
        breaks = ba_breaks_orig,
        labels = ba_breaks_labels
      ) +
      scale_y_continuous(
        trans = pseudo_log_trans(sigma = 1),
        breaks = c(0, 1, 2, 4, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = ba_x_limits, ylim = y_limits) +
      labs(
        x = expression("Pinyon Pine BA (m"^2*"/ha)"),
        y = if (is_first_col) "Expected Count" else NULL,
        title = trt_label
      ) +
      theme_bw(base_size = base_font_size) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = if (is_last_col) "right" else "none",
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        legend.key.width = unit(1.5, "cm")
      )
    
    plots_pinba[[trt]] <- p_pinba
    
    # Create ponderosa pine BA plot (flipped)
    p_pipoba <- ggplot(preds_pipoba_combined,
                       aes(x = x_original, y = .epred,
                           color = model_id, linetype = prop_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = model_id),
                  alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors, name = "Model") +
      scale_fill_manual(values = model_colors, guide = "none") +
      scale_linetype_manual(
        values = setNames(
          c("solid", "dashed", "dotted", "longdash")[1:length(prop_labels)],
          prop_labels
        ),
        name = "% Treated",
        drop = FALSE
      ) +
      scale_x_continuous(
        trans = pseudo_log_trans(sigma = 0.5, base = exp(1)),
        breaks = ba_breaks_orig,
        labels = ba_breaks_labels
      ) +
      scale_y_continuous(
        trans = pseudo_log_trans(sigma = 1),
        breaks = c(0, 1, 2, 4, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = ba_x_limits, ylim = y_limits) +
      labs(
        x = expression("Ponderosa Pine BA (m"^2*"/ha)"),
        y = if (is_first_col) "Expected Count" else NULL,
        title = trt_label
      ) +
      theme_bw(base_size = base_font_size) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = if (is_last_col) "right" else "none",
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        legend.key.width = unit(1.5, "cm")
      )
    
    plots_pipoba[[trt]] <- p_pipoba
  }
  
  # Combine predictions
  all_preds_df <- bind_rows(all_preds)
  
  # Combined plots
  combined_yrsince <- wrap_plots(plots_yrsince, nrow = 1, guides = "collect")
  combined_pinba <- wrap_plots(plots_pinba, nrow = 1, guides = "collect")
  combined_pipoba <- wrap_plots(plots_pipoba, nrow = 1, guides = "collect")
  
  # Combined BA plot with both pine species (2 rows)
  combined_ba <- wrap_plots(
    combined_pinba,
    combined_pipoba,
    nrow = 2,
    guides = "collect",
    axes = "collect"
  ) +
    plot_annotation(
      title = "Basal Area Interactions",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
    )
  
  return(list(
    predictions = all_preds_df,
    plots_yrsince = plots_yrsince,
    plots_pinba = plots_pinba,
    plots_pipoba = plots_pipoba,
    combined_plot_yrsince = combined_yrsince,
    combined_plot_pinba = combined_pinba,
    combined_plot_pipoba = combined_pipoba,
    combined_plot_ba = combined_ba
  ))
}


# Plot a single treatment's interaction effects across models

# Simplified interface for comparing 2 models on a single treatment
# Creates a figure with years-since and BA panels side by side or stacked
plot_treatment_comparison <- function(
    models,
    original_data_list,
    treatment,
    model_colors = NULL,
    plot_type = "flipped",
    interaction_var = "both",
    ...
) {
  
  if (plot_type == "flipped") {
    results <- plot_zip_interaction_effects_flipped_multimodel(
      models = models,
      original_data_list = original_data_list,
      treatments = treatment,
      model_colors = model_colors,
      ...
    )
  } else {
    results <- plot_zip_interaction_effects_multimodel(
      models = models,
      original_data_list = original_data_list,
      treatments = treatment,
      model_colors = model_colors,
      ...
    )
  }
  
  if (interaction_var == "yrsince") {
    return(results$plots_yrsince[[treatment]])
  } else if (interaction_var == "pinba") {
    return(results$plots_pinba[[treatment]])
  } else {
    # Both - combine vertically
    p <- results$plots_yrsince[[treatment]] / results$plots_pinba[[treatment]] +
      plot_layout(guides = "collect")
    return(p)
  }
}



