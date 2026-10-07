compare_zip_effects_multimodel <- function(
    models,  # Named list of models, e.g., list("Breeding" = m1, "Non-breeding" = m2)
    csv_output = "zip_effects_comparison.csv",
    plot_output = "zip_effects_comparison.png",
    effects_to_include = NULL,  # NULL = all effects, or character vector of effect names
    effect_labels = NULL,  # Named vector for relabeling, e.g., c("old_name" = "New Label")
    effect_groups = NULL,  # Named list of effect groupings for subplots, e.g., list("Controls" = c("year_scaled", ...), "Fire" = c(...))
    group_xlim = NULL,  # Named list of xlim for each group, e.g., list("Controls" = c(-1, 1), "Fire" = c(-2, 2))
    shared_xlim_groups = NULL,  # Character vector of group names that should share x-axis limits
    model_colors = NULL,  # Named vector of colors for models
    plot_width = 10,
    plot_height = 8,
    ci_level = 0.90,
    xlim = NULL,  # Global xlim (used when no grouping or as default)
    dodge_width = 0.6,
    point_size = 3,
    line_width = 1,
    show_component_prefix = TRUE,
    include_zi_effects = TRUE,
    base_font_size = 14,
    tag_style = "a)",  # Options: "a)", "a", "A)", "A", "(a)", "(A)", "1)", "1", "(1)"
    ncol = 1,  # Number of columns for subplot layout
    subplot_heights = NULL,  # Relative heights for subplots (NULL = auto based on n effects)
    legend_position = "top"  # Where to place the combined legend
) {
  
  require(brms)
  require(tidybayes)
  require(tidyverse)
  require(posterior)
  require(patchwork)
  
  # Validate input
  if (!is.list(models) || is.null(names(models))) {
    stop("'models' must be a named list, e.g., list('Full Year' = m1, 'Breeding' = m2)")
  }
  
  model_names <- names(models)
  message("Comparing ", length(models), " models: ", paste(model_names, collapse = ", "))
  
  # Helper function to extract effects from a single model
  extract_single_model <- function(model, model_name, ci_level) {
    
    draws <- as_draws_df(model)
    
    count_effects <- names(draws)[grepl("^b_(?!zi_)", names(draws), perl = TRUE)]
    zi_effects <- names(draws)[grepl("^b_zi_", names(draws))]
    
    count_effects <- count_effects[count_effects != "b_Intercept"]
    zi_effects <- zi_effects[zi_effects != "b_zi_Intercept"]
    
    all_effects <- c(count_effects, zi_effects)
    
    if (length(all_effects) == 0) {
      warning("No fixed effects found in model: ", model_name)
      return(NULL)
    }
    
    summary_stats <- map_dfr(all_effects, function(effect) {
      samples <- draws[[effect]]
      
      mean_val <- mean(samples)
      median_val <- median(samples)
      
      lower_ci <- (1 - ci_level) / 2
      upper_ci <- 1 - lower_ci
      ci_lower <- quantile(samples, lower_ci)
      ci_upper <- quantile(samples, upper_ci)
      
      p_positive <- mean(samples > 0)
      p_negative <- mean(samples < 0)
      
      component <- ifelse(grepl("^b_zi_", effect), "Zero-Inflation", "Count")
      effect_clean <- gsub("^b_zi_|^b_", "", effect)
      
      tibble(
        model_id = model_name,
        component = component,
        parameter = effect_clean,
        mean = mean_val,
        median = median_val,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        p_positive = p_positive,
        p_negative = p_negative
      )
    })
    
    return(summary_stats)
  }
  
  # Extract effects from all models
  message("Extracting fixed effects from all models...")
  all_results <- map_dfr(model_names, function(nm) {
    message("  Processing: ", nm)
    extract_single_model(models[[nm]], nm, ci_level)
  })
  
  if (nrow(all_results) == 0) {
    stop("No effects extracted from any model")
  }
  
  # Apply custom labels if provided
  if (!is.null(effect_labels)) {
    all_results <- all_results %>%
      mutate(parameter_label = ifelse(
        parameter %in% names(effect_labels),
        effect_labels[parameter],
        parameter
      ))
  } else {
    all_results <- all_results %>%
      mutate(parameter_label = parameter)
  }
  
  # Export to CSV
  message("Exporting results to: ", csv_output)
  write_csv(all_results, csv_output)
  
  # Prepare plot data
  message("Creating comparison dotplot...")
  
  plot_data <- all_results
  
  # Filter out ZI effects if requested
  if (!include_zi_effects) {
    plot_data <- plot_data %>%
      filter(component == "Count")
    message("Excluding zero-inflation effects from plot")
  }
  
  # Filter to specified effects
  if (!is.null(effects_to_include)) {
    plot_data <- plot_data %>%
      filter(parameter %in% effects_to_include)
    
    if (nrow(plot_data) == 0) {
      warning("No effects matched for plotting. Available effects: ", 
              paste(unique(all_results$parameter), collapse = ", "))
      plot_data <- all_results
    } else {
      message("Plotting ", length(unique(plot_data$parameter)), " specified effects")
      plot_data <- plot_data %>%
        mutate(parameter = factor(parameter, levels = effects_to_include))
    }
  }
  
  # Create plot labels
  plot_data <- plot_data %>%
    mutate(
      plot_label = if (show_component_prefix) {
        paste0("[", component, "] ", parameter_label)
      } else {
        parameter_label
      }
    )
  
  # Set model as factor
  plot_data <- plot_data %>%
    mutate(model_id = factor(model_id, levels = model_names))
  
  # Set up colors
  if (is.null(model_colors)) {
    default_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
    model_colors <- setNames(default_colors[1:length(model_names)], model_names)
  }
  
  # Parse tag style
  parse_tag_style <- function(style) {
    style <- tolower(style)
    if (grepl("^\\(", style)) {
      # Parentheses style: (a), (A), (1)
      prefix <- "("
      suffix <- ")"
      if (grepl("1", style)) {
        level <- "1"
      } else {
        level <- "a"
      }
    } else if (grepl("\\)$", style)) {
      # Suffix paren: a), A), 1)
      prefix <- ""
      suffix <- ")"
      if (grepl("1", style)) {
        level <- "1"
      } else {
        level <- "a"
      }
    } else {
      # No parentheses
      prefix <- ""
      suffix <- ""
      if (grepl("1", style)) {
        level <- "1"
      } else {
        level <- "a"
      }
    }
    list(level = level, prefix = prefix, suffix = suffix)
  }
  
  # Helper function to create a single subplot
  create_subplot <- function(data, group_name = NULL, xlim_vals = NULL, 
                             show_legend = FALSE, show_y_title = TRUE) {
    
    # Set effect order within this group
    if (!is.null(effects_to_include) && !is.null(effect_groups)) {
      # Use order from effect_groups
      group_effects <- effect_groups[[group_name]]
      label_order <- data %>%
        mutate(param_factor = factor(parameter, levels = group_effects)) %>%
        arrange(param_factor) %>%
        pull(plot_label) %>%
        unique()
    } else if (!is.null(effects_to_include)) {
      label_order <- data %>%
        arrange(parameter) %>%
        pull(plot_label) %>%
        unique()
    } else {
      label_order <- data %>%
        filter(model_id == model_names[1]) %>%
        arrange(component, median) %>%
        pull(plot_label) %>%
        unique()
    }
    
    data <- data %>%
      mutate(plot_label = factor(plot_label, levels = label_order))
    
    p <- ggplot(data, aes(x = median, y = plot_label, color = model_id)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      geom_point(
        position = position_dodge(width = dodge_width),
        size = point_size
      ) +
      geom_errorbarh(
        aes(xmin = ci_lower, xmax = ci_upper),
        position = position_dodge(width = dodge_width),
        height = 0,
        linewidth = line_width
      ) +
      scale_color_manual(
        values = model_colors,
        name = "Model"
      ) +
      labs(
        x = "Effect Size (standardized)",
        y = NULL
      ) +
      theme_bw(base_size = base_font_size) +
      theme(
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      )
    
    # Apply xlim
    if (!is.null(xlim_vals)) {
      p <- p + coord_cartesian(xlim = xlim_vals)
    }
    
    return(p)
  }
  
  # Determine if we're using grouped subplots or a single plot
  if (!is.null(effect_groups)) {
    
    message("Creating ", length(effect_groups), " grouped subplots...")
    
    # Validate effect_groups
    group_names <- names(effect_groups)
    if (is.null(group_names)) {
      stop("effect_groups must be a named list")
    }
    
    # Create list to store subplots
    subplot_list <- list()
    
    # Determine which groups share scales
    if (is.null(shared_xlim_groups)) {
      # By default, all groups with the same xlim share scales
      shared_xlim_groups <- group_names
    }
    
    # Calculate shared xlim for groups that should share
    shared_xlim_data <- plot_data %>%
      filter(parameter %in% unlist(effect_groups[shared_xlim_groups]))
    
    if (nrow(shared_xlim_data) > 0 && is.null(xlim)) {
      # Auto-calculate shared xlim with some padding
      data_range <- range(c(shared_xlim_data$ci_lower, shared_xlim_data$ci_upper), na.rm = TRUE)
      padding <- diff(data_range) * 0.1
      auto_shared_xlim <- c(data_range[1] - padding, data_range[2] + padding)
    } else {
      auto_shared_xlim <- xlim
    }
    
    # Create each subplot
    for (i in seq_along(effect_groups)) {
      group_name <- group_names[i]
      group_effects <- effect_groups[[group_name]]
      
      message("  Creating subplot: ", group_name, " (", length(group_effects), " effects)")
      
      # Filter data for this group
      group_data <- plot_data %>%
        filter(parameter %in% group_effects)
      
      if (nrow(group_data) == 0) {
        warning("No data for group: ", group_name)
        next
      }
      
      # Determine xlim for this group
      if (!is.null(group_xlim) && group_name %in% names(group_xlim)) {
        this_xlim <- group_xlim[[group_name]]
      } else if (group_name %in% shared_xlim_groups) {
        this_xlim <- auto_shared_xlim
      } else {
        this_xlim <- xlim
      }
      
      # Create subplot (no legend - we'll add a combined one)
      subplot_list[[group_name]] <- create_subplot(
        data = group_data,
        group_name = group_name,
        xlim_vals = this_xlim,
        show_legend = FALSE
      ) 
    }
    
    # Calculate subplot heights based on number of effects if not provided
    if (is.null(subplot_heights)) {
      subplot_heights <- sapply(effect_groups, length)
      # Normalize to reasonable range
      subplot_heights <- subplot_heights / max(subplot_heights)
      # Add minimum height
      subplot_heights <- pmax(subplot_heights, 0.3)
    }
    
    # Combine subplots with patchwork
    tag_info <- parse_tag_style(tag_style)
    
    # Create the combined plot
    if (ncol == 1) {
      # Vertical stack
      combined_plot <- wrap_plots(subplot_list, ncol = 1, heights = subplot_heights)
    } else {
      # Grid layout
      combined_plot <- wrap_plots(subplot_list, ncol = ncol)
    }
    
    # Add tags
    combined_plot <- combined_plot +
      plot_annotation(
        tag_levels = tag_info$level,
        tag_prefix = tag_info$prefix,
        tag_suffix = tag_info$suffix
      ) &
      theme(plot.tag = element_text(face = "bold", size = base_font_size + 2))
    
    # Add shared legend using patchwork's guide_area() approach
    # First, add legend back to the first subplot for extraction
    subplot_with_legend <- subplot_list[[1]] +
      theme(legend.position = legend_position)
    
    # Use patchwork's collect approach
    if (legend_position %in% c("top", "bottom")) {
      # Recreate subplots with legends, then collect
      subplot_list_with_legend <- subplot_list
      subplot_list_with_legend[[1]] <- subplot_with_legend
      
      if (ncol == 1) {
        combined_plot <- wrap_plots(subplot_list_with_legend, ncol = 1, heights = subplot_heights) +
          plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") +
          plot_annotation(
            tag_levels = tag_info$level,
            tag_prefix = tag_info$prefix,
            tag_suffix = tag_info$suffix
          ) &
          theme(
            plot.tag = element_text(face = "bold", size = base_font_size + 2),
            legend.position = legend_position
          )
      } else {
        combined_plot <- wrap_plots(subplot_list_with_legend, ncol = ncol) +
          plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") +
          plot_annotation(
            tag_levels = tag_info$level,
            tag_prefix = tag_info$prefix,
            tag_suffix = tag_info$suffix
          ) &
          theme(
            plot.tag = element_text(face = "bold", size = base_font_size + 2),
            legend.position = legend_position
          )
      }
    } else {
      # For right/left legend positions
      subplot_list_with_legend <- subplot_list
      subplot_list_with_legend[[1]] <- subplot_with_legend
      
      combined_plot <- wrap_plots(subplot_list_with_legend, ncol = ncol) +
        plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") +
        plot_annotation(
          tag_levels = tag_info$level,
          tag_prefix = tag_info$prefix,
          tag_suffix = tag_info$suffix
        ) &
        theme(
          plot.tag = element_text(face = "bold", size = base_font_size + 2),
          legend.position = legend_position
        )
    }
    
    p <- combined_plot
    
  } else {
    # Single plot (original behavior)
    message("Creating single comparison plot...")
    
    # Set plot order for effects
    if (!is.null(effects_to_include)) {
      label_order <- plot_data %>%
        arrange(parameter) %>%
        pull(plot_label) %>%
        unique()
      plot_data <- plot_data %>%
        mutate(plot_label = factor(plot_label, levels = label_order))
    } else {
      label_order <- plot_data %>%
        filter(model_id == model_names[1]) %>%
        arrange(component, median) %>%
        pull(plot_label) %>%
        unique()
      plot_data <- plot_data %>%
        mutate(plot_label = factor(plot_label, levels = label_order))
    }
    
    p <- ggplot(plot_data, aes(x = median, y = plot_label, color = model_id)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      geom_point(
        position = position_dodge(width = dodge_width),
        size = point_size
      ) +
      geom_errorbarh(
        aes(xmin = ci_lower, xmax = ci_upper),
        position = position_dodge(width = dodge_width),
        height = 0,
        linewidth = line_width
      ) +
      scale_color_manual(
        values = model_colors,
        name = "Model"
      ) +
      labs(
        x = "Effect Size (standardized)",
        y = NULL
      ) +
      theme_bw(base_size = base_font_size) +
      theme(
        legend.position = legend_position,
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      )
    
    if (!is.null(xlim)) {
      p <- p + coord_cartesian(xlim = xlim)
    }
  }
  
  # Save plot
  message("Saving plot to: ", plot_output)
  ggsave(plot_output, p, width = plot_width, height = plot_height, units = "in", dpi = 300)
  
  # Summary messages
  message("\nSummary:")
  message("  CSV contains ", nrow(all_results), " total effect estimates")
  message("  Plot contains ", nrow(plot_data), " effect estimates")
  for (nm in model_names) {
    n_effects <- sum(all_results$model_id == nm)
    message("  ", nm, ": ", n_effects, " effects")
  }
  if (!is.null(effect_groups)) {
    message("  Subplots created: ", length(effect_groups))
  }
  
  return(invisible(all_results))
}
