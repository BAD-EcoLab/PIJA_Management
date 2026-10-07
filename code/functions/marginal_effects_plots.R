## Marginal effects plots for ZIP model - Non-interacting fixed effects
## Uses tidybayes to extract posterior draws and plot expected counts
## Incorporates zero-inflation for overall expected values

library(tidyverse)
library(brms)
library(tidybayes)
library(ggdist)
library(modelr)
library(patchwork)
library(scales)  # For label_number()

plot_zip_marginal_effects <- function(
    model, #Fitted brms ZIP model
    original_data, #Original UNSCALED data (for back-transformation of axes)
    predictors = c("year_scaled", "day_of_year_scaled", "effort_hrs_scaled",
                   "num_observers_scaled", "cci_scaled", 
                   "elevation_30m_median_scaled", "propwater_scaled",
                   "pin_ba_log_scaled"), #Character vector of predictor names to plot (scaled versions)
    predictor_labels = NULL, #Named vector for nice axis labels, e.g., c("year_scaled" = "Year")
    n_points = 51, #Number of points along predictor range for predictions
    ci_width = 0.95, #Width of credible interval ribbon (default 0.95)
    n_draws = NULL #Number of posterior draws to use (NULL = all)
) {
  
  # Get the scaled data from the model
  model_data <- model$data
  
  # Identify all scaled predictors in the model data
  all_scaled_vars <- names(model_data)[grepl("_scaled$", names(model_data))]
  
  # Binary treatment indicators - set to 0 (no treatment) for baseline
  binary_vars <- c("prescribed_fire_bin", "mastication_bin", "thinning_bin",
                   "harvest_bin", "is_mobile")
  
  # Create a baseline data frame with means for scaled vars, 0 for binary
  baseline <- model_data %>%
    summarise(across(all_of(all_scaled_vars), mean, na.rm = TRUE))
  
  # Add binary treatment variables set to 0 (no treatment baseline)
  for (bv in binary_vars) {
    baseline[[bv]] <- 0
  }
  
  message("Baseline: scaled predictors at mean, treatment indicators = 0 (no treatment)")
  
  # Calculate scaling parameters from original data for back-transformation
  scaling_params <- list()
  for (pred in predictors) {
    orig_name <- gsub("_scaled$", "", pred)
    if (orig_name %in% names(original_data)) {
      scaling_params[[pred]] <- list(
        mean = mean(original_data[[orig_name]], na.rm = TRUE),
        sd = sd(original_data[[orig_name]], na.rm = TRUE),
        orig_name = orig_name
      )
    }
  }
  
  # Identify predictors with quadratic effects in the model
  # These need their squared terms updated when varying the linear term
  quadratic_predictors <- c("Iday_of_year_scaledE2", "Ielevation_30m_median_scaledE2", 
                            "pin_ba_log_scaledE2", "pipo_ba_log_scaledE2")
  
  # Default labels if not provided
  if (is.null(predictor_labels)) {
    predictor_labels <- setNames(
      c("Year", "Day of Year", "Effort (hrs)", "# Observers", 
        "CCI", "Elevation (m)", "Prop. Water", 
        "Pinyon BA (log)", "Ponderosa Pine BA (log)"),
      predictors
    )
  }
  
  # Store predictions for each predictor
  all_preds <- list()
  plots <- list()
  
  for (pred in predictors) {
    message("Processing: ", pred)
    
    # Get the range of this predictor in the model data
    pred_range <- range(model_data[[pred]], na.rm = TRUE)
    
    # Create prediction grid: vary this predictor, hold others at baseline
    pred_grid <- baseline %>%
      slice(rep(1, n_points)) %>%
      mutate(!!pred := seq(pred_range[1], pred_range[2], length.out = n_points))
    
    # For predictors with quadratic effects, update the squared term
    # This handles models with I(x^2) or pre-computed squared columns
    if (pred %in% quadratic_predictors) {
      sq_name <- paste0(pred, "_sq")
      # Check if squared column exists in model data (pre-computed)
      if (sq_name %in% names(model_data)) {
        pred_grid <- pred_grid %>%
          mutate(!!sq_name := .data[[pred]]^2)
      }
      # Note: If model uses I(x^2) in formula, brms handles it automatically
    }
    
    # Add ecoregion for random effect structure
    pred_grid$ecoregion <- "new"
    
    # Extract posterior predictions - OVERALL expected value (includes zero-inflation)
    # By NOT specifying dpar, we get E[Y] = (1 - zi) * lambda
    preds <- pred_grid %>%
      add_epred_draws(
        model, 
        re_formula = NA,  # Exclude random effects (population-level prediction)
        allow_new_levels = TRUE,
        ndraws = n_draws
        # No dpar argument = overall expected value incorporating ZI
      )
    
    # Summarize across draws
    pred_summary <- preds %>%
      group_by(across(all_of(pred))) %>%
      median_qi(.epred, .width = ci_width) %>%
      ungroup() %>%
      mutate(predictor = pred)
    
    # Back-transform predictor to original scale
    if (pred %in% names(scaling_params)) {
      sp <- scaling_params[[pred]]
      pred_summary <- pred_summary %>%
        mutate(x_original = .data[[pred]] * sp$sd + sp$mean)
    } else {
      pred_summary <- pred_summary %>%
        mutate(x_original = .data[[pred]])
    }
    
    all_preds[[pred]] <- pred_summary
    
    # Get axis label
    x_label <- ifelse(pred %in% names(predictor_labels), 
                      predictor_labels[pred], 
                      gsub("_scaled$", "", pred))
    
    # Determine if this is a left-column plot (for y-axis label logic)
    pred_index <- which(predictors == pred)
    is_left_column <- (pred_index %% 2) == 1  # Odd indices are left column
    
    # Create individual plot with original-scale x-axis
    p <- ggplot(pred_summary, aes(x = x_original, y = .epred)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper), 
                  fill = "steelblue", alpha = 0.3) +
      geom_line(color = "steelblue", linewidth = 1) +
      scale_x_continuous(labels = scales::label_number(accuracy = 1)) +  # Integer labels
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +  # Start at 0
      labs(
        x = x_label,
        y = if (is_left_column) "Expected Abundance" else NULL  # Only left column gets y label
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      )
    
    plots[[pred]] <- p
  }
  
  # Combine all predictions
  all_preds_df <- bind_rows(all_preds)
  
  # Create combined plot using patchwork
  message("Creating combined plot...")
  combined <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = "Marginal Effects: Non-interacting Predictors",
      subtitle = paste0("Median with ", ci_width * 100, "% CI | Other predictors at mean, no treatment")
    )
  
  return(list(
    predictions = all_preds_df,
    plots = plots,
    combined_plot = combined
  ))
}


#' Alternative plotting using stat_lineribbon (closer to example code style)
#' Shows nested credible intervals
plot_marginal_lineribbon <- function(
    model,
    original_data,
    predictor,
    predictor_label = NULL,
    n_points = 51,
    ci_widths = c(0.5, 0.8, 0.95),
    n_draws = NULL,
    color = "steelblue"
) {
  
  model_data <- model$data
  
  # Get baseline (same logic as above)
  all_scaled_vars <- names(model_data)[grepl("_scaled$", names(model_data))]
  binary_vars <- c("prescribed_fire_bin", "mastication_bin", "thinning_bin",
                   "harvest_bin", "is_mobile")# "clearcut_bin")
  
  baseline <- model_data %>%
    summarise(across(all_of(all_scaled_vars), mean, na.rm = TRUE))
  
  # Add binary treatment variables set to 0
  for (bv in binary_vars) {
    baseline[[bv]] <- 0
  }
  
  pred_range <- range(model_data[[predictor]], na.rm = TRUE)
  
  pred_grid <- baseline %>%
    slice(rep(1, n_points)) %>%
    mutate(!!predictor := seq(pred_range[1], pred_range[2], length.out = n_points))
  
  # For predictors with quadratic effects, update the squared term
  quadratic_predictors <- c("year_scaled", "elevation_30m_median_scaled", "pin_ba_log_scaled")
  if (predictor %in% quadratic_predictors) {
    sq_name <- paste0(predictor, "_sq")
    if (sq_name %in% names(model_data)) {
      pred_grid <- pred_grid %>%
        mutate(!!sq_name := .data[[predictor]]^2)
    }
  }
  
  pred_grid$ecoregion <- "new"
  
  # Get draws (not summarized) - overall expected value
  preds <- pred_grid %>%
    add_epred_draws(
      model,
      re_formula = NA,
      allow_new_levels = TRUE,
      ndraws = n_draws
      # No dpar = overall expected value
    )
  
  # Back-transform predictor (individual scaling)
  orig_name <- gsub("_scaled$", "", predictor)
  if (orig_name %in% names(original_data)) {
    orig_mean <- mean(original_data[[orig_name]], na.rm = TRUE)
    orig_sd <- sd(original_data[[orig_name]], na.rm = TRUE)
    preds <- preds %>%
      mutate(x_original = .data[[predictor]] * orig_sd + orig_mean)
  } else {
    preds <- preds %>%
      mutate(x_original = .data[[predictor]])
  }
  
  # Axis label
  if (is.null(predictor_label)) {
    predictor_label <- gsub("_scaled$", "", predictor)
  }
  
  # Plot with stat_lineribbon
  p <- ggplot(preds, aes(x = x_original, y = .epred)) +
    stat_lineribbon(.width = ci_widths, fill = color, alpha = 0.4) +
    scale_fill_brewer(palette = "Blues") +
    scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = predictor_label,
      y = "Expected Abundance"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  return(p)
}


#' Plot interaction effects for treatment variables in ZIP model
#'
#' Creates two separate figures:
#' 1. Time x Treatment: Effects by years since treatment
#' 2. BA x Treatment: Effects by Pinyon and Ponderosa basal area (2 rows)
#'
#' @param model Fitted brms ZIP model
#' @param original_data Original UNSCALED data (for back-transformation)
#' @param treatments Character vector of treatment names (without _scaled suffix)
#' @param treatment_labels Named vector for nice treatment labels
#' @param yrsince_levels Years since treatment to plot (default: c(1, 5, 10))
#' @param pinba_levels Absolute basal area values for pinyon pine (will be log-transformed with 0.1 offset)
#' @param pipoba_levels Absolute basal area values for ponderosa pine (will be log-transformed with 0.1 offset)
#' @param n_points Number of points along treatment proportion range
#' @param ci_width Width of credible interval ribbon
#' @param n_draws Number of posterior draws to use (NULL = all)
#'
#' @return A list with: predictions, plots, combined_plot_yrsince, combined_plot_ba
plot_zip_interaction_effects <- function(
    model,
    original_data,
    treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
    treatment_labels = NULL,
    yrsince_levels = c(1, 5, 10),
    pinba_levels = c(1, 5, 10),
    pipoba_levels = c(1, 5, 10),
    n_points = 51,
    ci_width = 0.9,
    n_draws = NULL
) {
  
  model_data <- model$data
  
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
  
  # Get all scaled variables and binary indicators
  all_scaled_vars <- names(model_data)[grepl("_scaled$", names(model_data))]
  binary_vars <- c("prescribed_fire_bin", "mastication_bin", "thinning_bin",
                   "harvest_bin", "is_mobile")
  
  # Calculate scaling parameters for back-transformation and forward-transformation
  # Using individual scaling for each variable
  scaling_params <- list()
  
  for (trt in treatments) {
    # Treatment proportion
    if (trt %in% names(original_data)) {
      scaling_params[[paste0(trt, "_scaled")]] <- list(
        mean = mean(original_data[[trt]], na.rm = TRUE),
        sd = sd(original_data[[trt]], na.rm = TRUE)
      )
    }
    # Years since
    yrsince_var <- paste0(trt, "_yrsince")
    if (yrsince_var %in% names(original_data)) {
      scaling_params[[paste0(yrsince_var, "_scaled")]] <- list(
        mean = mean(original_data[[yrsince_var]], na.rm = TRUE),
        sd = sd(original_data[[yrsince_var]], na.rm = TRUE)
      )
    }
  }
  
  # Pin BA scaling
  if ("pin_ba_log" %in% names(original_data)) {
    scaling_params[["pin_ba_log_scaled"]] <- list(
      mean = mean(original_data$pin_ba_log, na.rm = TRUE),
      sd = sd(original_data$pin_ba_log, na.rm = TRUE)
    )
    # Convert absolute BA values to log scale (with 0.1 offset to avoid log(0))
    pinba_values_orig <- log(pinba_levels + 0.1)
    pinba_labels <- as.character(pinba_levels)
  }
  
  # Pipo BA scaling
  if ("pipo_ba_log" %in% names(original_data)) {
    scaling_params[["pipo_ba_log_scaled"]] <- list(
      mean = mean(original_data$pipo_ba_log, na.rm = TRUE),
      sd = sd(original_data$pipo_ba_log, na.rm = TRUE)
    )
    # Convert absolute BA values to log scale (with 0.1 offset to avoid log(0))
    pipoba_values_orig <- log(pipoba_levels + 0.1)
    pipoba_labels <- as.character(pipoba_levels)
  }
  
  # Create baseline with all scaled vars at mean, all treatments OFF
  baseline <- model_data %>%
    summarise(across(all_of(all_scaled_vars), mean, na.rm = TRUE))
  for (bv in binary_vars) {
    baseline[[bv]] <- 0
  }
  
  message("Creating interaction effect plots...")
  
  # Storage for all predictions and plots
  all_preds <- list()
  plots_yrsince <- list()
  plots_pinba <- list()
  plots_pipoba <- list()
  
  for (trt in treatments) {
    message("Processing treatment: ", trt)
    
    trt_scaled <- paste0(trt, "_scaled")
    trt_bin <- paste0(trt, "_bin")
    yrsince_scaled <- paste0(trt, "_yrsince_scaled")
    
    # Get range of treatment proportion (in scaled units)
    trt_range <- range(model_data[[trt_scaled]], na.rm = TRUE)
    
    # Get label
    trt_label <- ifelse(trt %in% names(treatment_labels),
                        treatment_labels[trt], trt)
    
    # Determine column position for axis label logic
    trt_index <- which(treatments == trt)
    is_first_col <- trt_index == 1
    is_last_col <- trt_index == length(treatments)
    
    # =========================================================================
    # PLOT 1: Varying years since treatment
    # =========================================================================
    
    sp_yrsince <- scaling_params[[yrsince_scaled]]
    yrsince_scaled_values <- (yrsince_levels - sp_yrsince$mean) / sp_yrsince$sd
    
    preds_yrsince <- map_dfr(seq_along(yrsince_levels), function(i) {
      pred_grid <- baseline %>%
        slice(rep(1, n_points)) %>%
        mutate(
          !!trt_scaled := seq(trt_range[1], trt_range[2], length.out = n_points),
          !!trt_bin := 1,
          !!yrsince_scaled := yrsince_scaled_values[i]
        )
      pred_grid$ecoregion <- "new"
      
      preds <- pred_grid %>%
        add_epred_draws(
          model,
          re_formula = NA,
          allow_new_levels = TRUE,
          ndraws = n_draws
        ) %>%
        group_by(across(all_of(trt_scaled))) %>%
        median_qi(.epred, .width = ci_width) %>%
        ungroup() %>%
        mutate(
          yrsince_raw = yrsince_levels[i],
          yrsince_label = paste0(yrsince_levels[i], " yr")
        )
      return(preds)
    })
    
    sp_trt <- scaling_params[[trt_scaled]]
    preds_yrsince <- preds_yrsince %>%
      mutate(
        x_original = .data[[trt_scaled]] * sp_trt$sd + sp_trt$mean,
        treatment = trt,
        # Order labels by numeric year value, not alphabetically
        yrsince_label = factor(yrsince_label, levels = paste0(sort(yrsince_levels), " yr"))
      )
    
    all_preds[[paste0(trt, "_yrsince")]] <- preds_yrsince
    
    # Create years_since plot (with treatment label as title for column header)
    p_yrs <- ggplot(preds_yrsince, aes(x = x_original, y = .epred, 
                                       color = yrsince_label, fill = yrsince_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_viridis_d(option = "B", end = 0.8, name = "Years Since") +
      scale_fill_viridis_d(option = "B", end = 0.8, name = "Years Since") +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 1),
        breaks = c(0, 1, 3, 6, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 10)) +
      labs(
        x = "Proportion Treated",
        y = if (is_first_col) "Relative Abundance" else NULL,
        title = trt_label
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = if (is_last_col) "right" else "none",
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold")
      )
    
    plots_yrsince[[trt]] <- p_yrs
    
    # =========================================================================
    # PLOT 2: Varying pinyon pine basal area (TOP ROW of BA figure)
    # =========================================================================
    
    sp_pinba <- scaling_params[["pin_ba_log_scaled"]]
    pinba_scaled_values <- (pinba_values_orig - sp_pinba$mean) / sp_pinba$sd
    
    preds_pinba <- map_dfr(seq_along(pinba_levels), function(i) {
      pred_grid <- baseline %>%
        slice(rep(1, n_points)) %>%
        mutate(
          !!trt_scaled := seq(trt_range[1], trt_range[2], length.out = n_points),
          !!trt_bin := 1,
          pin_ba_log_scaled = pinba_scaled_values[i]
        )
      pred_grid$ecoregion <- "new"
      
      preds <- pred_grid %>%
        add_epred_draws(
          model,
          re_formula = NA,
          allow_new_levels = TRUE,
          ndraws = n_draws
        ) %>%
        group_by(across(all_of(trt_scaled))) %>%
        median_qi(.epred, .width = ci_width) %>%
        ungroup() %>%
        mutate(pinba_label = pinba_labels[i])
      return(preds)
    })
    
    preds_pinba <- preds_pinba %>%
      mutate(
        x_original = .data[[trt_scaled]] * sp_trt$sd + sp_trt$mean,
        treatment = trt,
        pinba_label = factor(pinba_label, levels = pinba_labels)
      )
    
    all_preds[[paste0(trt, "_pinba")]] <- preds_pinba
    
    # Create pin_ba plot (TOP ROW of BA figure - no title, column labels from yrsince)
    p_pinba <- ggplot(preds_pinba, aes(x = x_original, y = .epred,
                                       color = pinba_label, fill = pinba_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_viridis_d(option = "B", end = 0.8, name = "Pinyon BA") +
      scale_fill_viridis_d(option = "B", end = 0.8, name = "Pinyon BA") +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 1),
        breaks = c(0, 1, 3, 6, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 10)) +
      labs(
        x = NULL,  # No x-axis title for top row
        y = if (is_first_col) "Relative Abundance" else NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = if (is_last_col) "right" else "none",
        axis.text.x = element_blank(),  # No x-axis tick labels for top row
        axis.ticks.x = element_blank()
      )
    
    plots_pinba[[trt]] <- p_pinba
    
    # =========================================================================
    # PLOT 3: Varying ponderosa pine basal area (BOTTOM ROW of BA figure)
    # =========================================================================
    
    sp_pipoba <- scaling_params[["pipo_ba_log_scaled"]]
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
        add_epred_draws(
          model,
          re_formula = NA,
          allow_new_levels = TRUE,
          ndraws = n_draws
        ) %>%
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
        pipoba_label = factor(pipoba_label, levels = pipoba_labels)
      )
    
    all_preds[[paste0(trt, "_pipoba")]] <- preds_pipoba
    
    # Create pipo_ba plot (BOTTOM ROW - has x-axis labels and title)
    p_pipoba <- ggplot(preds_pipoba, aes(x = x_original, y = .epred,
                                         color = pipoba_label, fill = pipoba_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_viridis_d(option = "B", end = 0.8, name = "Ponderosa BA") +
      scale_fill_viridis_d(option = "B", end = 0.8, name = "Ponderosa BA") +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 1),
        breaks = c(0, 1, 3, 6, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 10)) +
      labs(
        x = "Proportion Treated",
        y = if (is_first_col) "Relative Abundance" else NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = if (is_last_col) "right" else "none"
      )
    
    plots_pipoba[[trt]] <- p_pipoba
  }
  
  # Combine all predictions
  all_preds_df <- bind_rows(all_preds)
  
  # ==========================================================================
  # CREATE COMBINED PLOTS
  # ==========================================================================
  
  message("Creating combined plots...")
  
  # FIGURE 1: Years since treatment (single row)
  # Treatment labels appear as titles on each subplot
  combined_yrsince <- wrap_plots(plots_yrsince, nrow = 1)
  
  # FIGURE 2: BA interactions (two rows: pinyon on top, ponderosa on bottom)
  # No subplot titles - column labels inherited from yrsince row when combined
  row_pinba <- wrap_plots(plots_pinba, nrow = 1)
  row_pipoba <- wrap_plots(plots_pipoba, nrow = 1)
  
  combined_ba <- wrap_plots(
    row_pinba,
    row_pipoba,
    nrow = 2,
    guides = "collect",
    axes = "collect"
  )
  
  return(list(
    predictions = all_preds_df,
    plots_yrsince = plots_yrsince,
    plots_pinba = plots_pinba,
    plots_pipoba = plots_pipoba,
    combined_plot_yrsince = combined_yrsince,
    combined_plot_ba = combined_ba
  ))
}


#' Plot interaction effects with TIME/BA on x-axis and proportion treated as groups
#'
#' This is the "flipped" version of plot_zip_interaction_effects:
#' - Row 1: Years since treatment (1-10) on x-axis, lines for different proportion treated levels
#' - Row 2: Pinyon BA (0 to max) on x-axis, lines for different proportion treated levels
#' - Row 3: Ponderosa BA (0 to max) on x-axis, lines for different proportion treated levels
#'
#' @param model Fitted brms ZIP model
#' @param original_data Original UNSCALED data (for back-transformation)
#' @param treatments Character vector of treatment names (without _scaled suffix)
#' @param treatment_labels Named vector for nice treatment labels
#' @param prop_treated_levels Proportion treated levels to show (default: c(0.25, 0.50, 0.75))
#' @param yrsince_range Range of years since treatment for x-axis (default: c(1, 10))
#' @param n_points Number of points along x-axis range
#' @param ci_width Width of credible interval ribbon
#' @param n_draws Number of posterior draws to use (NULL = all)
#'
#' @return A list with: predictions, plots, combined plots
plot_zip_interaction_effects_flipped <- function(
    model,
    original_data,
    treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
    treatment_labels = NULL,
    prop_treated_levels = c(0.25, 0.50, 0.75),
    yrsince_range = c(1, 10),
    n_points = 51,
    ci_width = 0.9,
    n_draws = NULL
) {
  
  model_data <- model$data
  
  # Default treatment labels
  if (is.null(treatment_labels)) {
    treatment_labels <- c(
      "prescribed_fire" = "Prescribed Fire",
      "mastication" = "Mastication", 
      "thinning" = "Thinning",
      "harvest" = "Harvest"
      # "clearcut" = "Clearcut"
    )
  }
  
  # Labels for proportion treated levels
  prop_labels <- paste0(prop_treated_levels * 100, "%")
  
  # Get all scaled variables and binary indicators
  all_scaled_vars <- names(model_data)[grepl("_scaled$", names(model_data))]
  binary_vars <- c("prescribed_fire_bin", "mastication_bin", "thinning_bin",
                   "harvest_bin", "is_mobile")
  
  # Calculate scaling parameters for back-transformation and forward-transformation
  scaling_params <- list()
  max_prop_observed <- list()  # Store max observed proportion for each treatment
  
  for (trt in treatments) {
    # Treatment proportion
    if (trt %in% names(original_data)) {
      scaling_params[[paste0(trt, "_scaled")]] <- list(
        mean = mean(original_data[[trt]], na.rm = TRUE),
        sd = sd(original_data[[trt]], na.rm = TRUE)
      )
      # Store max observed proportion for this treatment
      max_prop_observed[[trt]] <- max(original_data[[trt]], na.rm = TRUE)
    }
    # Years since
    yrsince_var <- paste0(trt, "_yrsince")
    if (yrsince_var %in% names(original_data)) {
      scaling_params[[paste0(yrsince_var, "_scaled")]] <- list(
        mean = mean(original_data[[yrsince_var]], na.rm = TRUE),
        sd = sd(original_data[[yrsince_var]], na.rm = TRUE)
      )
    }
  }
  
  # Pin BA scaling and range
  # Note: pin_ba_log is already log-transformed, so we use min/max of log values
  if ("pin_ba_log" %in% names(original_data)) {
    scaling_params[["pin_ba_log_scaled"]] <- list(
      mean = mean(original_data$pin_ba_log, na.rm = TRUE),
      sd = sd(original_data$pin_ba_log, na.rm = TRUE)
    )
    # Use actual data range (log scale) - includes values where original BA < 1
    pinba_range_orig <- c(min(original_data$pin_ba_log, na.rm = TRUE), 
                          max(original_data$pin_ba_log, na.rm = TRUE))
  }
  
  # Pipo BA scaling and range
  if ("pipo_ba_log" %in% names(original_data)) {
    scaling_params[["pipo_ba_log_scaled"]] <- list(
      mean = mean(original_data$pipo_ba_log, na.rm = TRUE),
      sd = sd(original_data$pipo_ba_log, na.rm = TRUE)
    )
    # Use actual data range (log scale)
    pipoba_range_orig <- c(min(original_data$pipo_ba_log, na.rm = TRUE),
                           max(original_data$pipo_ba_log, na.rm = TRUE))
  }
  
  # X-axis limits on ORIGINAL BA scale (for pseudo-log transformation)
  ba_xlim_orig <- c(0, 100)
  
  # Break points for BA on original scale
  # Pseudo-log is linear near zero, logarithmic for larger values
  ba_breaks_orig <- c(0, 1, 5, 25, 100)
  ba_breaks_labels <- c("0", "1", "5", "25", "100")
  
  # Create baseline with all scaled vars at mean, all treatments OFF
  baseline <- model_data %>%
    summarise(across(all_of(all_scaled_vars), mean, na.rm = TRUE))
  for (bv in binary_vars) {
    baseline[[bv]] <- 0
  }
  
  message("Creating flipped interaction effect plots...")
  message("X-axis: years since / basal area | Color: proportion treated levels")
  
  # Storage for all predictions and plots
  all_preds <- list()
  plots_yrsince <- list()
  plots_pinba <- list()
  plots_pinba_titled <- list()  # Version with titles for BA-only figure
  plots_pipoba <- list()
  
  for (trt in treatments) {
    message("Processing treatment: ", trt)
    
    trt_scaled <- paste0(trt, "_scaled")
    trt_bin <- paste0(trt, "_bin")
    yrsince_scaled <- paste0(trt, "_yrsince_scaled")
    
    # Get scaling params for this treatment
    sp_trt <- scaling_params[[trt_scaled]]
    sp_yrsince <- scaling_params[[yrsince_scaled]]
    
    # Filter prop_treated_levels to only those within observed range for this treatment
    max_prop <- max_prop_observed[[trt]]
    valid_prop_idx <- which(prop_treated_levels <= max_prop)
    
    if (length(valid_prop_idx) == 0) {
      message("  Skipping ", trt, ": no prop_treated_levels within observed range (max = ", 
              round(max_prop, 2), ")")
      next
    }
    
    # Filter to valid levels
    trt_prop_levels <- prop_treated_levels[valid_prop_idx]
    trt_prop_labels <- paste0(trt_prop_levels * 100, "%")
    trt_prop_scaled <- (trt_prop_levels - sp_trt$mean) / sp_trt$sd
    
    message("  Using prop levels: ", paste(trt_prop_labels, collapse = ", "), 
            " (max observed = ", round(max_prop, 2), ")")
    
    # Get label
    trt_label <- ifelse(trt %in% names(treatment_labels),
                        treatment_labels[trt], trt)
    
    # Determine column position for axis label logic
    trt_index <- which(treatments == trt)
    is_first_col <- trt_index == 1
    is_last_col <- trt_index == length(treatments)
    
    # =========================================================================
    # PLOT 1: Years since treatment on x-axis
    # =========================================================================
    
    # Convert yrsince range to scaled values
    yrsince_range_scaled <- (yrsince_range - sp_yrsince$mean) / sp_yrsince$sd
    
    preds_yrsince <- map_dfr(seq_along(trt_prop_levels), function(i) {
      pred_grid <- baseline %>%
        slice(rep(1, n_points)) %>%
        mutate(
          !!yrsince_scaled := seq(yrsince_range_scaled[1], yrsince_range_scaled[2], length.out = n_points),
          !!trt_scaled := trt_prop_scaled[i],
          !!trt_bin := 1
        )
      pred_grid$ecoregion <- "new"
      
      preds <- pred_grid %>%
        add_epred_draws(
          model,
          re_formula = NA,
          allow_new_levels = TRUE,
          ndraws = n_draws
        ) %>%
        group_by(across(all_of(yrsince_scaled))) %>%
        median_qi(.epred, .width = ci_width) %>%
        ungroup() %>%
        mutate(
          prop_raw = trt_prop_levels[i],
          prop_label = trt_prop_labels[i]
        )
      return(preds)
    })
    
    # Back-transform years since to original scale
    preds_yrsince <- preds_yrsince %>%
      mutate(
        x_original = .data[[yrsince_scaled]] * sp_yrsince$sd + sp_yrsince$mean,
        treatment = trt,
        # Use FULL prop_labels as factor levels for consistent colors across plots
        prop_label = factor(prop_label, levels = prop_labels)
      )
    
    all_preds[[paste0(trt, "_yrsince_flipped")]] <- preds_yrsince
    
    # Create years_since plot (with treatment label as title for first row)
    p_yrs <- ggplot(preds_yrsince, aes(x = x_original, y = .epred, 
                                       color = prop_label, fill = prop_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_viridis_d(option = "B", end = 0.8, name = "% Treated", drop = FALSE) +
      scale_fill_viridis_d(option = "B", end = 0.8, name = "% Treated", drop = FALSE) +
      scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 1),
        breaks = c(0, 1, 2, 4, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(ylim = c(-0.1, 4)) +
      labs(
        x = "Years Since Treatment",
        # y = if (is_first_col) "Relative Abundance" else NULL,
        y = "Relative Abundance"
        # title = trt_label  # Treatment label at top of column
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank()
        # legend.position = if (is_last_col) "right" else "none",
        # plot.title = element_text(hjust = 0.5, size = 11, face = "bold")
      )
    
    plots_yrsince[[trt]] <- p_yrs
    
    # =========================================================================
    # PLOT 2: Pinyon pine basal area on x-axis
    # =========================================================================
    
    sp_pinba <- scaling_params[["pin_ba_log_scaled"]]
    pinba_range_scaled <- (pinba_range_orig - sp_pinba$mean) / sp_pinba$sd
    
    preds_pinba <- map_dfr(seq_along(trt_prop_levels), function(i) {
      pred_grid <- baseline %>%
        slice(rep(1, n_points)) %>%
        mutate(
          pin_ba_log_scaled = seq(pinba_range_scaled[1], pinba_range_scaled[2], length.out = n_points),
          !!trt_scaled := trt_prop_scaled[i],
          !!trt_bin := 1
        )
      pred_grid$ecoregion <- "new"
      
      preds <- pred_grid %>%
        add_epred_draws(
          model,
          re_formula = NA,
          allow_new_levels = TRUE,
          ndraws = n_draws
        ) %>%
        group_by(pin_ba_log_scaled) %>%
        median_qi(.epred, .width = ci_width) %>%
        ungroup() %>%
        mutate(
          prop_raw = trt_prop_levels[i],
          prop_label = trt_prop_labels[i]
        )
      return(preds)
    })
    
    # Back-transform pin BA to original scale (exp of log values)
    preds_pinba <- preds_pinba %>%
      mutate(
        # Back-transform: first unscale to log values, then exp() to original BA
        x_original = exp(pin_ba_log_scaled * sp_pinba$sd + sp_pinba$mean),
        treatment = trt,
        # Use FULL prop_labels as factor levels for consistent colors across plots
        prop_label = factor(prop_label, levels = prop_labels)
      )
    
    all_preds[[paste0(trt, "_pinba_flipped")]] <- preds_pinba
    
    # Create pinyon BA plot (no title - for combined_all where yrsince has column labels)
    p_pinba <- ggplot(preds_pinba, aes(x = x_original, y = .epred,
                                        color = prop_label, fill = prop_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_viridis_d(option = "B", end = 0.8, name = "% Treated", drop = FALSE) +
      scale_fill_viridis_d(option = "B", end = 0.8, name = "% Treated", drop = FALSE) +
      scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.5, base = exp(1)),
                         breaks = ba_breaks_orig, labels = ba_breaks_labels) +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 1),
        breaks = c(0, 1, 2, 4, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = ba_xlim_orig, 
                      ylim = c(-0.1, 6)) +
      labs(
        x = expression("Pinyon Pine BA (m"^2*"/ha)"),
        y = "Relative Abundance"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank()
      )
    
    plots_pinba[[trt]] <- p_pinba
    
    # Create titled version for BA-only figure (where pinba is top row)
    p_pinba_titled <- p_pinba +
      labs(title = trt_label) +
      theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
    
    plots_pinba_titled[[trt]] <- p_pinba_titled
    
    # =========================================================================
    # PLOT 3: Ponderosa pine basal area on x-axis
    # =========================================================================
    
    sp_pipoba <- scaling_params[["pipo_ba_log_scaled"]]
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
        add_epred_draws(
          model,
          re_formula = NA,
          allow_new_levels = TRUE,
          ndraws = n_draws
        ) %>%
        group_by(pipo_ba_log_scaled) %>%
        median_qi(.epred, .width = ci_width) %>%
        ungroup() %>%
        mutate(
          prop_raw = trt_prop_levels[i],
          prop_label = trt_prop_labels[i]
        )
      return(preds)
    })
    
    # Back-transform pipo BA to original scale (exp of log values)
    preds_pipoba <- preds_pipoba %>%
      mutate(
        # Back-transform: first unscale to log values, then exp() to original BA
        x_original = exp(pipo_ba_log_scaled * sp_pipoba$sd + sp_pipoba$mean),
        treatment = trt,
        # Use FULL prop_labels as factor levels for consistent colors across plots
        prop_label = factor(prop_label, levels = prop_labels)
      )
    
    all_preds[[paste0(trt, "_pipoba_flipped")]] <- preds_pipoba
    
    # Create ponderosa BA plot (no title - column labels from row above)
    p_pipoba <- ggplot(preds_pipoba, aes(x = x_original, y = .epred,
                                          color = prop_label, fill = prop_label)) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_viridis_d(option = "B", end = 0.8, name = "% Treated", drop = FALSE) +
      scale_fill_viridis_d(option = "B", end = 0.8, name = "% Treated", drop = FALSE) +
      scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.5, base = exp(1)),
                         breaks = ba_breaks_orig, labels = ba_breaks_labels) +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 1),
        breaks = c(0, 1, 2, 4, 10),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = ba_xlim_orig, 
                      ylim = c(-0.1, 6)) +
      labs(
        x = expression("Ponderosa Pine BA (m"^2*"/ha)"),
        y = 'Relative Abundance'
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank()
      )
    
    plots_pipoba[[trt]] <- p_pipoba
  }
  
  # Combine all predictions
  all_preds_df <- bind_rows(all_preds)
  
  # ==========================================================================
  # CREATE COMBINED PLOTS
  # ==========================================================================
  
  message("Creating combined plots...")
  
  # FIGURE 1: Years since treatment on x-axis (single row)
  # Treatment labels appear as titles on each subplot
  # combined_yrsince <- wrap_plots(plots_yrsince, nrow = 1)
  
  # Figire 1: rx fire, mastication, thinning; single column
  combined_yrsince <- 
    plots_yrsince$prescribed_fire + 
    annotate('text', x = 10, y = 3.8, hjust = 1,
             label = "Prescribed Fire",
             fontface = "bold") +
    plots_yrsince$mastication +
      annotate('text', x = 10, y = 3.8, hjust = 1,
               label = "Mastication",
               fontface = "bold") +
    plots_yrsince$thinning +
      annotate('text', x = 10, y = 3.5, hjust = 1,
               label = "Thinning",
               fontface = "bold") +
    plot_layout(guides = "collect", axes = 'collect',
                nrow = 3) +
    plot_annotation(tag_levels = 'a', tag_suffix = ')') 

  
  # FIGURE 2: BA on x-axis (two rows: pinyon on top, ponderosa on bottom)
  # Use titled version of pinba for column labels since this is standalone
  # row_pinba_titled <- wrap_plots(plots_pinba_titled, nrow = 1)
  row_pipoba <- wrap_plots(plots_pipoba, nrow = 1)
  # 
  # combined_ba <- row_pinba_titled / row_pipoba
  
  combined_ba <- 
    plots_pinba$prescribed_fire + 
    annotate('text', x = 100, y = 6, hjust = 1,
             label = "Prescribed Fire",
             fontface = "bold") +
    plots_pinba$mastication +
    annotate('text', x = 100, y = 6, hjust = 1,
             label = "Mastication",
             fontface = "bold") +
    plots_pinba$thinning +
    annotate('text', x = 100, y = 6, hjust = 1,
             label = "Thinning",
             fontface = "bold") +
    plots_pinba$harvest + theme(legend.position = 'none') +
    annotate('text', x = 100, y = 6, hjust = 1,
             label = "Harvest",
             fontface = "bold") +
    plots_pipoba$prescribed_fire +
    plots_pipoba$mastication +
    plots_pipoba$thinning +
    plots_pipoba$harvest + theme(legend.position = 'none') +
    plot_layout(guides = "collect", axes = 'collect',
                nrow = 2)
  
  # FIGURE 3: All three rows combined
  # Treatment labels on top row (yrsince), BA rows underneath (no titles)
  row_yrs <- wrap_plots(plots_yrsince, nrow = 1)
  row_pinba <- wrap_plots(plots_pinba, nrow = 1)
  combined_all <- row_yrs / row_pinba / row_pipoba
  
  return(list(
    predictions = all_preds_df,
    plots_yrsince = plots_yrsince,
    plots_pinba = plots_pinba,
    plots_pinba_titled = plots_pinba_titled,
    plots_pipoba = plots_pipoba,
    combined_plot_yrsince = combined_yrsince,
    combined_plot_ba = combined_ba,
    combined_plot_all = combined_all
  ))
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# # Load your fitted model and original data
# fit <- readRDS("path/to/pinjay_zip.rds")
# original_data <- arrow::read_parquet("path/to/original_data.parquet")

# -----------------------------------------------------------------------------
# NON-INTERACTING EFFECTS
# -----------------------------------------------------------------------------

# # Basic usage - all 8 non-interacting predictors
# results_main <- plot_zip_marginal_effects(
#   model = fit,
#   original_data = original_data,
#   n_points = 51,
#   ci_width = 0.95,
#   n_draws = 500  # Use subset of draws for speed; NULL for all
# )
# 
# # View and save
# results_main$combined_plot
# ggsave("marginal_effects_main.png", results_main$combined_plot,
#        width = 10, height = 12, dpi = 300)

# -----------------------------------------------------------------------------
# TREATMENT INTERACTION EFFECTS
# -----------------------------------------------------------------------------

# # Plot treatment interactions with years_since, pin_ba, and pipo_ba
# results_interactions <- plot_zip_interaction_effects(
#   model = fit,
#   original_data = original_data,
#   treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
#   yrsince_levels = c(1, 5, 10),  # Years since treatment to show
#   pinba_levels = c(1, 5, 10),  # Absolute pinyon BA values (will be log-transformed)
#   pipoba_levels = c(1, 5, 10),  # Absolute ponderosa BA values (will be log-transformed)
#   n_points = 51,
#   ci_width = 0.95,
#   n_draws = 500
# )
# 
# # View and save the years-since figure (single row)
# results_interactions$combined_plot_yrsince
# ggsave("treatment_yrsince.png", results_interactions$combined_plot_yrsince,
#        width = 14, height = 4, dpi = 300)
#
# # View and save the BA figure (two rows: pinyon on top, ponderosa on bottom)
# results_interactions$combined_plot_ba
# ggsave("treatment_ba.png", results_interactions$combined_plot_ba,
#        width = 14, height = 8, dpi = 300)
# 
# # Access individual plots
# results_interactions$plots_yrsince$prescribed_fire
# results_interactions$plots_pinba$thinning
# results_interactions$plots_pipoba$mastication

# -----------------------------------------------------------------------------
# FLIPPED INTERACTION EFFECTS (time/BA on x-axis, % treated as groups)
# -----------------------------------------------------------------------------

# # Plot with years-since and BA on x-axis, different lines for % treated
# results_flipped <- plot_zip_interaction_effects_flipped(
#   model = fit,
#   original_data = original_data,
#   treatments = c("prescribed_fire", "mastication", "thinning", "harvest"),
#   prop_treated_levels = c(0.25, 0.50, 0.75),  # 25%, 50%, 75% treated
#   yrsince_range = c(1, 10),  # Years since treatment range for x-axis
#   n_points = 51,
#   ci_width = 0.95,
#   n_draws = 500
# )
# 
# # View and save the years-since figure (time on x-axis)
# results_flipped$combined_plot_yrsince
# ggsave("treatment_over_time.png", results_flipped$combined_plot_yrsince,
#        width = 14, height = 4, dpi = 300)
#
# # View and save the BA figure (BA on x-axis, two rows)
# results_flipped$combined_plot_ba
# ggsave("treatment_across_ba.png", results_flipped$combined_plot_ba,
#        width = 14, height = 8, dpi = 300)
#
# # View all three rows combined
# results_flipped$combined_plot_all
# ggsave("treatment_interactions_flipped.png", results_flipped$combined_plot_all,
#        width = 14, height = 10, dpi = 300)
# 
# # Custom proportion treated levels (e.g., 10%, 30%, 50%, 70%)
# results_custom <- plot_zip_interaction_effects_flipped(
#   model = fit,
#   original_data = original_data,
#   prop_treated_levels = c(0.10, 0.30, 0.50, 0.70),
#   n_draws = 500
# )

# -----------------------------------------------------------------------------
# CUSTOM LABELS
# -----------------------------------------------------------------------------

# # Custom treatment labels
# custom_trt_labels <- c(
#   "prescribed_fire" = "Rx Fire",
#   "mastication" = "Mastication",
#   "thinning" = "Thinning",
#   "harvest" = "Harvest",
#   "clearcut" = "Clearcut"
# )
# 
# results_interactions <- plot_zip_interaction_effects(
#   model = fit,
#   original_data = original_data,
#   treatment_labels = custom_trt_labels,
#   n_draws = 500
# )


# =============================================================================
# QUICK START: If you don't have original data, you can approximate
# =============================================================================

# # If original data is unavailable, the model$data has standardized values
# # The back-transformation will just show standardized units
# # OR you can manually provide the scaling parameters:
# 
# # Example: create a fake "original_data" with known ranges
# # (You would replace these with actual values from your data processing)
# approx_original <- tibble(
#   year = seq(2010, 2023, length.out = 100),
#   day_of_year = seq(1, 365, length.out = 100),
#   effort_hrs = seq(0.1, 12, length.out = 100),
#   num_observers = seq(1, 20, length.out = 100),
#   cci = seq(0, 1, length.out = 100),
#   elevation_30m_median = seq(1000, 3500, length.out = 100),
#   propwater = seq(0, 0.5, length.out = 100),
#   pin_ba_log = seq(-2, 4, length.out = 100),
#   pipo_ba_log = seq(-2, 4, length.out = 100)
# )


# =============================================================================
# CUSTOMIZATION OPTIONS
# =============================================================================

# # Custom axis labels
# custom_labels <- c(
#   "year_scaled" = "Year",
#   "day_of_year_scaled" = "Day of Year", 
#   "effort_hrs_scaled" = "Survey Effort (hours)",
#   "num_observers_scaled" = "Number of Observers",
#   "cci_scaled" = "Cloud Cover Index",
#   "elevation_30m_median_scaled" = "Elevation (m)",
#   "propwater_scaled" = "Proportion Water",
#   "pin_ba_log_scaled" = "Pinyon Basal Area (log)"
# )
# 
# results <- plot_zip_marginal_effects(
#   model = fit,
#   original_data = original_data,
#   predictor_labels = custom_labels,
#   n_draws = 1000
# )


# =============================================================================
# LINERIBBON STYLE (like your example code)
# =============================================================================

# # Individual plots with nested credible intervals
# p_elev <- plot_marginal_lineribbon(
#   model = fit,
#   original_data = original_data,
#   predictor = "elevation_30m_median_scaled",
#   predictor_label = "Elevation (m)",
#   color = "steelblue"
# )
# 
# p_year <- plot_marginal_lineribbon(
#   model = fit,
#   original_data = original_data,
#   predictor = "year_scaled",
#   predictor_label = "Year",
#   color = "darkgreen"
# )
# 
# # Combine with patchwork
# p_elev + p_year


# =============================================================================
# OUTPUT STRUCTURE
# =============================================================================

# The plot_zip_marginal_effects function returns a list with:
# - predictions: tibble with columns:
#     - [predictor]: scaled predictor value
#     - .epred: median expected abundance
#     - .lower, .upper: CI bounds
#     - predictor: name of the predictor
#     - x_original: back-transformed predictor value
# - plots: named list of individual ggplots
# - combined_plot: patchwork combined plot

# The plot_zip_interaction_effects function returns a list with:
# - predictions: tibble with all predictions
# - plots_yrsince: named list of years-since interaction plots (x = prop treated)
# - plots_pinba: named list of pinyon BA interaction plots (x = prop treated)
# - plots_pipoba: named list of ponderosa BA interaction plots (x = prop treated)
# - combined_plot_yrsince: patchwork combined plot for time interactions
# - combined_plot_ba: patchwork combined plot for BA interactions (2 rows)

# The plot_zip_interaction_effects_flipped function returns a list with:
# - predictions: tibble with all predictions
# - plots_yrsince: named list (x = years since, color = % treated)
# - plots_pinba: named list (x = pinyon BA, color = % treated)
# - plots_pipoba: named list (x = ponderosa BA, color = % treated)
# - combined_plot_yrsince: single row showing time effects
# - combined_plot_ba: two rows showing BA effects (pinyon, ponderosa)
# - combined_plot_all: three rows (years, pinyon BA, ponderosa BA)