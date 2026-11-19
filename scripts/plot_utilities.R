# All the plot functions used for the rumination LNM project
# Xiao Chen 260627
# chenxiaophd@gmail.com

############# spatial corr bar plot ########################
plot_correlation_bar <- function(
    df, out_file,
    annotation_levels = c("margulies2016_fcgradient01_fsLR_32k",
                          "margulies2016_fcgradient02_fsLR_32k",
                          "margulies2016_fcgradient03_fsLR_32k"),
    annotation_labels = c("Gradient 1", "Gradient 2", "Gradient 3"),
    fill_colors = c("#FF894F", "#FFCB61", "#77BEF0"),
    width_mm = 55, height_mm = 100,
    y_label_fontsize = 9,
    x_tick_fontsize = 9,
    y_tick_fontsize = 8
) {
  # Axis title and tick label font sizes are in points.
  
  library(ggplot2)
  library(dplyr)
  
  df <- df %>%
    mutate(
      sig_label = ifelse(p < 0.05, "**", ""),
      annotation = factor(
        annotation,
        levels = annotation_levels
      )
    )
  
  p <- ggplot(
    df,
    aes(
      x = annotation,
      y = r,
      fill = annotation
    )
  ) +
    
    geom_col(
      width = 0.55,
      colour = NA
    ) +
    
    geom_text(
      aes(label = sig_label),
      vjust = ifelse(df$r >= 0, -0.12, 0.9),
      size = 6,
      fontface = "bold",
      family = "Arial"
    ) +
    
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      colour = "black"
    ) +
    
    scale_fill_manual(
      values = fill_colors
    ) +
    
    scale_x_discrete(
      labels = annotation_labels
    ) +
    
    labs(
      x = NULL,
      y = expression("Spatial corr (" * italic(r) * ")")
    ) +
    
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 9) +
    theme(
      text = element_text(family = "Arial",colour = "black"),
      axis.line = element_line(linewidth = 0.5,colour = "black"),
      axis.ticks = element_line(linewidth = 0.5,colour = "black"),
      axis.text.x = element_text(angle = 50,
        hjust = 1,
        vjust = 1,
        size = x_tick_fontsize,
        colour = "black",
        family = "Arial"
      ),
      
      axis.text.y = element_text(
        size = y_tick_fontsize,
        colour = "black"
      ),
      
      axis.title.y = element_text(
        size = y_label_fontsize,
        colour = "black"
      ),
      
      legend.position = "none",
      
      plot.margin = margin(5, 8, 5, 5)
    )
  
  ggsave(
    out_file,
    p,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    dpi = 600
  )
  
  return(p)
}

######################## word cloud plot ########################
plot_word_cloud <- function(
    df,
    out_file = NULL,
    label_col = "Terms",
    corr_col = "corr",
    seed = 45,
    vertical_probability = 0.35,
    first_n_horizontal = 3,
    low_color = "#90CAF9",
    high_color = "#0D47A1",
    size_range = c(6, 22),
    width = 35,
    height = 40,
    units = "cm",
    dpi = 300,
    top_n = NULL,
    size_multiplier = 1
) {
  library(ggplot2)
  library(dplyr)
  library(ggwordcloud)

  required_cols <- c(label_col, corr_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (
    length(vertical_probability) != 1 ||
    is.na(vertical_probability) ||
    vertical_probability < 0 ||
    vertical_probability > 1
  ) {
    stop("vertical_probability must be a single number between 0 and 1.")
  }

  if (
    length(first_n_horizontal) != 1 ||
    is.na(first_n_horizontal) ||
    first_n_horizontal < 0 ||
    first_n_horizontal %% 1 != 0
  ) {
    stop("first_n_horizontal must be a non-negative integer.")
  }

  if (nrow(df) == 0L || !is.numeric(df[[corr_col]]) ||
      any(!is.finite(df[[corr_col]]))) {
    stop("Correlation values must be numeric, finite, and non-empty.")
  }
  if (length(size_multiplier) != 1L || !is.finite(size_multiplier) ||
      size_multiplier <= 0) {
    stop("size_multiplier must be a positive finite number.")
  }
  # Keep the original size/color mapping when selecting fewer words.
  size_limits <- range(pmax(df[[corr_col]], 0))
  color_limits <- range(df[[corr_col]])
  if (!is.null(top_n)) {
    if (length(top_n) != 1L || !is.finite(top_n) ||
        top_n < 1 || top_n %% 1 != 0) {
      stop("top_n must be a positive integer.")
    }
    df <- df %>% arrange(desc(.data[[corr_col]])) %>% head(top_n)
  }

  set.seed(seed)
  plot_df <- df %>%
    mutate(
      angle = 90 * sample(
        c(0, 1),
        n(),
        replace = TRUE,
        prob = c(1 - vertical_probability, vertical_probability)
      ),
      corr_wc = pmax(.data[[corr_col]], 0)
#      corr_wc = sqrt(pmax(.data[[corr_col]], 0))
    )

  n_horizontal <- min(first_n_horizontal, nrow(plot_df))
  if (n_horizontal > 0) {
    plot_df$angle[seq_len(n_horizontal)] <- 0
  }

  p <- ggplot(
    plot_df,
    aes(
      label = .data[[label_col]],
      size = corr_wc,
      angle = angle,
      color = .data[[corr_col]]
    )
  ) +
    geom_text_wordcloud_area(
      shape = "circle",
      area_corr = TRUE,
      eccentricity = 1.5,
      rm_outside = FALSE,
      grid_size = 1,
      seed = seed,
      show.legend = FALSE
    ) +
    theme_minimal(base_size = 25) +
#    scale_radius(range = size_range, limits = c(0, NA)) +
    scale_radius(range = size_range * size_multiplier, limits = size_limits) +
    scale_color_gradient(
      name = "Corr",
      low = low_color,
      high = high_color,
      limits = color_limits
    ) +
    guides(size = "none")

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(
      filename = out_file,
      plot = p,
      width = width,
      height = height,
      units = units,
      dpi = dpi,
      bg = "white"
    )
  }

  return(p)
}

####################### binned LNM plot ######################
plot_gradient_binned_from_summary <- function(
    df_bins,
    summary_df,
    out_file = NULL,
    point_color = "#91ADC8",
    mean_color = "black",
    error = "sd",
    ylab = expression("LNM (" * italic(t) * ")"),
    xlab = "Gradient percentile"
){
  
  library(ggplot2)
  
  if (error == "sd") {
    summary_df$error <- summary_df$sd
  } else if (error == "sem") {
    summary_df$error <- summary_df$sem
  } else {
    stop("error must be 'sd' or 'sem'")
  }
  
  df_bins$bin <- factor(df_bins$bin, levels = sort(unique(df_bins$bin)))
  summary_df$bin <- factor(summary_df$bin, levels = levels(df_bins$bin))
  
  if (!"label" %in% colnames(summary_df)) {
    n_bins <- length(levels(df_bins$bin))
    summary_df$label <- paste0(
      seq(0, 100 - 100/n_bins, by = 100/n_bins),
      "-",
      seq(100/n_bins, 100, by = 100/n_bins),
      "%"
    )
  }
  
  p <- ggplot(df_bins, aes(x = bin, y = value)) +
    geom_hline(
      yintercept = 0,
      linetype = "22",      # dashed
      colour = "grey70",
      linewidth = 0.5
    )+
    geom_jitter(
      width = 0.18,
      size = 2,
      alpha = 0.5,
      shape = 21,
      fill = point_color,
      colour = "black",
      stroke = 0.1
    )+
    geom_errorbar(
      data = summary_df,
      aes(
        x = bin,
        ymin = mean - error,
        ymax = mean + error
      ),
      inherit.aes = FALSE,
      width = 0.13,
      linewidth = 0.7,
      colour = mean_color
    ) +
    geom_point(
      data = summary_df,
      aes(x = bin, y = mean),
      inherit.aes = FALSE,
      colour = mean_color,
      size = 3
    ) +
    scale_x_discrete(labels = summary_df$label) +
    labs(
      x = xlab,
      y = ylab
    ) +
    theme_classic(base_family = "Arial") +
    theme(
      axis.text = element_text(size = 16, colour = "black"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.title = element_text(size = 16) #face = "bold"
    )
  if (!is.null(out_file)) {
    ggsave(
      filename = out_file,
      plot = p,
      width = 6,
      height = 4,
      dpi = 600
    )
  }
  
  return(p)
}

##################### bar plot of neurotransmitter ###################
plot_nt_spin_bar <- function(
    df,
    out_file = NULL,
    top_n = NULL,
    sort_by = c("r_obs", "p_spin"),
    bar_color = "#77BEF0",
    width_mm = 170,
    x_label_fontsize = 16,
    x_tick_fontsize = 14,
    y_tick_fontsize = 14,
    height_mm = NULL
){
  # Axis title and tick label font sizes are in points.
  # height_mm is in millimetres; NULL uses max(100, number of plotted rows * 6).
  
  library(ggplot2)
  library(dplyr)
  
  sort_by <- match.arg(sort_by)
  
  df <- df %>%
    mutate(
      nt_label = sub("_.*$", "", label),
      
      sig_label = case_when(
        p_spin < 0.05 & p_fdr < 0.05 ~ "**",
        p_spin < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  if (!is.null(top_n)) {
    
    if (sort_by == "r_obs") {
      df <- df %>%
        arrange(desc(abs(r_obs))) %>%
        slice_head(n = top_n)
    }
    
    if (sort_by == "p_spin") {
      df <- df %>%
        arrange(p_spin) %>%
        slice_head(n = top_n)
    }
  }
  
  df <- df %>%
    arrange(r_obs) %>%
    mutate(
      nt_label = factor(nt_label, levels = nt_label)
    )
  
  offset <- 0.02 * diff(range(df$r_obs, na.rm = TRUE))
  
  p <- ggplot(
    df,
    aes(
      x = r_obs,
      y = nt_label
    )
  ) +
    
    geom_vline(
      xintercept = 0,
      colour = "grey70",
      linewidth = 0.4
    ) +
    
    geom_col(
      width = 0.72,
      colour = NA,
      fill = bar_color
    ) +
    
    geom_text(
      aes(
        x = ifelse(
          r_obs >= 0,
          r_obs + offset,
          r_obs - offset
        ),
        label = sig_label
      ),
      family = "Arial",
      fontface = "bold",
      size = 8,
      hjust = ifelse(df$r_obs >= 0, 0, 1),
      vjust = 0.75,
    ) +
    
    # scale_fill_manual(
    #   values = c(
    #     "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
    #     "#76B7B2", "#EDC948", "#B07AA1", "#FF9DA7",
    #     "#9C755F", "#BAB0AC", "#86BCB6", "#D37295",
    #     "#8CD17D", "#499894", "#D4A6C8", "#79706E",
    #     "#FABFD2", "#B6992D", "#A0CBE8", "#FFBE7D"
    #   )
    # )+
    
    labs(
      x = expression("Spatial correlation (" * italic(r) * ")"),
      y = NULL
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_classic(base_family = "Arial") +
    
    theme(
      text = element_text(
        size = 12,
        colour = "black"
      ),
      
      axis.title.x = element_text(
        size = x_label_fontsize
      ),
      
      axis.text.x = element_text(
        size = x_tick_fontsize,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        size = y_tick_fontsize,
        colour = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.45
      ),
      
      axis.ticks = element_line(
        linewidth = 0.45
      ),
      
      axis.ticks.length = unit(2, "mm"),
      
      legend.position = "none"
    )
  
  if (!is.null(out_file)) {
    ggsave(
      out_file,
      plot = p,
      width = width_mm,
      height = if (is.null(height_mm)) max(100, nrow(df) * 6) else height_mm,
      units = "mm",
      dpi = 600
    )
  }
  
  return(p)
}

###################### dot line plot ###############################
plot_two_group_paired <- function(
    df,
    value_col = "metric",
    condition_col = "condition",
    group_col = "group",
    subject_col = NULL,
    ylab = "ALFF",
    error = c("sd", "sem"),
    sig_df = NULL,
    out_file = NULL,
    width = 3.8,
    height = 3.2,
    dpi = 600,
    jitter_width = 0.1,
    jitter_seed = 123,
    xlab = NULL,
    x_tick_label_size = 12,
    y_tick_label_size = 12,
    x_axis_label_size = 12,
    y_axis_label_size = 12,
    group_label_size = 4 * 72.27 / 25.4
) {
  
  library(ggplot2)
  library(dplyr)
  
  error <- match.arg(error)
  df <- df %>%
    mutate(
      condition = factor(.data[[condition_col]]),
      group = factor(.data[[group_col]])
    )
  
  groups <- levels(df$group)
  conditions <- levels(df$condition)
  
  if (length(groups) != 2) stop("This function expects exactly 2 groups.")
  if (length(conditions) != 2) stop("This function expects exactly 2 conditions.")
  if (is.null(subject_col)) {
    n_cond <- length(conditions)
    df <- df %>%
      arrange(group, condition) %>%
      group_by(group) %>%
      mutate(
        subject_id = rep(seq_len(n() / n_cond), times = n_cond)
      ) %>%
      ungroup()
  } else {
    df <- df %>%
      mutate(subject_id = .data[[subject_col]])
  }
  df <- df %>%
    mutate(
      x_pos = case_when(
        group == groups[1] & condition == conditions[1] ~ 1,
        group == groups[1] & condition == conditions[2] ~ 2,
        group == groups[2] & condition == conditions[1] ~ 3.2,
        group == groups[2] & condition == conditions[2] ~ 4.2
      )
    )
  
  set.seed(jitter_seed)
  
  df <- df %>%
    mutate(
      x_jitter = x_pos + runif(n(), min = -jitter_width, max = jitter_width)
    )
  
  summary_df <- df %>%
    group_by(x_pos, group, condition) %>%
    summarise(
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      sd_value = sd(.data[[value_col]], na.rm = TRUE),
      n = sum(!is.na(.data[[value_col]])),
      .groups = "drop"
    ) %>%
    mutate(
      error_value = ifelse(
        error == "sd",
        sd_value,
        sd_value / sqrt(n)
      ),
      ymin = mean_value - error_value,
      ymax = mean_value + error_value
    )
  
  p <- ggplot(
    df,
    aes(x = x_jitter, y = .data[[value_col]])
  ) +
    geom_line(
      aes(group = interaction(group, subject_id)),
      colour = "grey78",
      linewidth = 0.28,
      alpha = 0.75
    ) +
    geom_point(
      aes(fill = condition),
      shape = 21,
      size = 3.5,
      stroke = 0.25,
      colour = "black",
      alpha = 0.85
    ) +
    geom_errorbar(
      data = summary_df,
      aes(
        x = x_pos,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      width = 0.10,
      linewidth = 0.45,
      colour = "black"
    ) +
    geom_segment(
      data = summary_df,
      aes(
        x = x_pos - 0.14,
        xend = x_pos + 0.14,
        y = mean_value,
        yend = mean_value
      ),
      inherit.aes = FALSE,
      linewidth = 0.65,
      colour = "black"
    )
  
  # -------------------------
  # Add significance annotation
  # -------------------------
  
  if (!is.null(sig_df) && nrow(sig_df) > 0) {
    p <- p +
      geom_segment(
        data = sig_df,
        aes(
          x = .data[["xmin"]],
          xend = .data[["xmax"]],
          y = .data[["y.position"]],
          yend = .data[["y.position"]]
        ),
        inherit.aes = FALSE,
        linewidth = 0.5,
        colour = "black"
      ) +
      geom_segment(
        data = sig_df,
        aes(
          x = .data[["xmin"]],
          xend = .data[["xmin"]],
          y = .data[["y.position"]],
          yend = .data[["y.position"]] - .data[["tip.length"]]
        ),
        inherit.aes = FALSE,
        linewidth = 0.5,
        colour = "black"
      ) +
      geom_segment(
        data = sig_df,
        aes(
          x = .data[["xmax"]],
          xend = .data[["xmax"]],
          y = .data[["y.position"]],
          yend = .data[["y.position"]] - .data[["tip.length"]]
        ),
        inherit.aes = FALSE,
        linewidth = 0.5,
        colour = "black"
      ) +
      geom_text(
        data = sig_df,
        aes(
          x = (.data[["xmin"]] + .data[["xmax"]]) / 2,
          y = .data[["y.position"]] + .data[["text.offset"]],
          label = .data[["label"]]
        ),
        inherit.aes = FALSE,
        size = 5,
        family = "Arial",
        colour = "black"
      )
  }
  
  p <- p +
    scale_x_continuous(
      breaks = c(1, 2, 3.2, 4.2),
      labels = c("Rum.", "Dis.", "Rum.", "Dis.")
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    annotate(
      "text",
      x = 1.5,
      y = -Inf,
      label = "MDD",
      vjust = 4,
      family = "Arial",
      size = group_label_size / (72.27 / 25.4)
    ) +
    annotate(
      "text",
      x = 3.7,
      y = -Inf,
      label = "HC",
      vjust = 4,
      family = "Arial",
      size = group_label_size / (72.27 / 25.4)
    ) +
    labs(
      x = xlab,
      y = ylab
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 7) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      axis.text.y = element_text(size = y_tick_label_size, colour = "black"),
      axis.text.x = element_text(
        size = x_tick_label_size,
        colour = "black",
        angle = 30,
        hjust = 1,
        vjust = 1.2
      ),
      axis.title.x = element_text(
        size = x_axis_label_size,
        colour = "black",
        margin = margin(t = 32)
      ),
      axis.title.y = element_text(
        size = y_axis_label_size,
        colour = "black",
        margin = margin(r = 4)
      ),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      legend.position = "none",
      plot.margin = margin(12, 4, 18, 4)
    )
  if (!is.null(out_file)) {
    ggsave(
      filename = out_file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  return(p)
}

#################### dot line plot, three groups ##################
plot_three_group_paired <- function(
    df,
    value_col = "metric",
    condition_col = "condition",
    group_col = "group",
    subject_col = NULL,
    ylab = "ALFF",
    error = c("sd", "sem"),
    sig_df = NULL,
    out_file = NULL,
    width = 4.5,
    height = 3.8,
    dpi = 600,
    jitter_width = 0.1,
    jitter_seed = 123
) {
  
  library(ggplot2)
  library(dplyr)
  
  error <- match.arg(error)
  
  df <- df %>%
    mutate(
      condition = factor(.data[[condition_col]]),
      group = factor(.data[[group_col]])
    )
  
  groups <- levels(df$group)
  conditions <- levels(df$condition)
  
  if (length(groups) != 3) stop("This function expects exactly 3 groups.")
  if (length(conditions) != 2) stop("This function expects exactly 2 conditions.")
  
  if (is.null(subject_col)) {
    n_cond <- length(conditions)
    
    df <- df %>%
      arrange(group, condition) %>%
      group_by(group) %>%
      mutate(
        subject_id = rep(seq_len(n() / n_cond), times = n_cond)
      ) %>%
      ungroup()
  } else {
    df <- df %>%
      mutate(subject_id = .data[[subject_col]])
  }
  
  # x positions:
  # group1: 1, 2
  # group2: 3.2, 4.2
  # group3: 5.4, 6.4
  
  x_map <- expand.grid(
    group = groups,
    condition = conditions,
    stringsAsFactors = FALSE
  ) %>%
    arrange(group, condition) %>%
    mutate(
      group_id = as.numeric(factor(group, levels = groups)),
      cond_id = as.numeric(factor(condition, levels = conditions)),
      x_pos = (group_id - 1) * 2.2 + cond_id
    )
  
  df <- df %>%
    left_join(
      x_map,
      by = c("group", "condition")
    )
  
  set.seed(jitter_seed)
  
  df <- df %>%
    mutate(
      x_jitter = x_pos + runif(
        n(),
        min = -jitter_width,
        max = jitter_width
      )
    )
  
  summary_df <- df %>%
    group_by(x_pos, group, condition) %>%
    summarise(
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      sd_value = sd(.data[[value_col]], na.rm = TRUE),
      n = sum(!is.na(.data[[value_col]])),
      .groups = "drop"
    ) %>%
    mutate(
      error_value = ifelse(
        error == "sd",
        sd_value,
        sd_value / sqrt(n)
      ),
      ymin = mean_value - error_value,
      ymax = mean_value + error_value
    )
  
  p <- ggplot(
    df,
    aes(
      x = x_jitter,
      y = .data[[value_col]]
    )
  ) +
    geom_line(
      aes(group = interaction(group, subject_id)),
      colour = "grey78",
      linewidth = 0.28,
      alpha = 0.75
    ) +
    geom_point(
      aes(fill = condition),
      shape = 21,
      size = 3.5,
      stroke = 0.25,
      colour = "black",
      alpha = 0.85
    ) +
    geom_errorbar(
      data = summary_df,
      aes(
        x = x_pos,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      width = 0.10,
      linewidth = 0.45,
      colour = "black"
    ) +
    geom_segment(
      data = summary_df,
      aes(
        x = x_pos - 0.14,
        xend = x_pos + 0.14,
        y = mean_value,
        yend = mean_value
      ),
      inherit.aes = FALSE,
      linewidth = 0.65,
      colour = "black"
    )
  
  if (!is.null(sig_df) && nrow(sig_df) > 0) {
    
    p <- p +
      geom_segment(
        data = sig_df,
        aes(
          x = .data[["xmin"]],
          xend = .data[["xmax"]],
          y = .data[["y.position"]],
          yend = .data[["y.position"]]
        ),
        inherit.aes = FALSE,
        linewidth = 0.5,
        colour = "black"
      ) +
      geom_segment(
        data = sig_df,
        aes(
          x = .data[["xmin"]],
          xend = .data[["xmin"]],
          y = .data[["y.position"]],
          yend = .data[["y.position"]] - .data[["tip.length"]]
        ),
        inherit.aes = FALSE,
        linewidth = 0.5,
        colour = "black"
      ) +
      geom_segment(
        data = sig_df,
        aes(
          x = .data[["xmax"]],
          xend = .data[["xmax"]],
          y = .data[["y.position"]],
          yend = .data[["y.position"]] - .data[["tip.length"]]
        ),
        inherit.aes = FALSE,
        linewidth = 0.5,
        colour = "black"
      ) +
      geom_text(
        data = sig_df,
        aes(
          x = (.data[["xmin"]] + .data[["xmax"]]) / 2,
          y = .data[["y.position"]] + .data[["text.offset"]],
          label = .data[["label"]]
        ),
        inherit.aes = FALSE,
        size = 5,
        family = "Arial",
        colour = "black"
      )
  }
  
  x_breaks <- x_map$x_pos
  x_labels <- rep(c("Rum.", "Dis."), times = length(groups))
  
  group_label_df <- x_map %>%
    group_by(group) %>%
    summarise(
      x = mean(x_pos),
      .groups = "drop"
    )
  
  p <- p +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    labs(
      x = NULL,
      y = ylab
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 7) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(
        size = 12,
        colour = "black",
        angle = 30,
        hjust = 1,
        vjust = 1.2
      ),
      axis.title.y = element_text(
        size = 12,
        colour = "black",
        margin = margin(r = 4)
      ),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      legend.position = "none",
      plot.margin = margin(12, 4, 22, 4)
    )
  
  for (i in seq_len(nrow(group_label_df))) {
    p <- p +
      annotate(
        "text",
        x = group_label_df$x[i],
        y = -Inf,
        label = group_label_df$group[i],
        vjust = 4,
        family = "Arial",
        size = 4
      )
  }
  
  if (!is.null(out_file)) {
    ggsave(
      filename = out_file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  return(p)
}

################### scatter plot ##########################
# Draw each plotmath line separately so unequal line lengths stay left-aligned.
element_grob.element_corr_stats <- function(element, label = "", ...) {
  class(element) <- setdiff(class(element), "element_corr_stats")
  rows <- lapply(as.list(label), function(line) {
    ggplot2::element_grob(element, label = as.expression(line), ...)
  })
  # gtable is already a required dependency of ggplot2.
  layout <- gtable::gtable(
    widths = grid::unit(1, "null"),
    heights = grid::unit(rep(element$size * 1.6, length(rows)), "pt")
  )
  gtable::gtable_add_grob(
    layout, grobs = rows, t = seq_along(rows), l = 1, clip = "off"
  )
}

plot_corr_txt <- function(
    df,
    out_file,
    title = "",
    xlab = "FC",
    ylab = "Affect Rating",
    p_adjust_factor = 6,
    x_axis_label_size = 30,
    y_axis_label_size = 30,
    x_tick_label_size = 24,
    y_tick_label_size = 24,
    stats_label_size = 18
) {
  
  library(ggplot2)
  
  df$Groups <- factor(
    df$Groups,
    levels = c(1, 2),
    labels = c("MDD", "HC")
  )

  if (length(p_adjust_factor) != 1 ||
      !is.finite(p_adjust_factor) ||
      p_adjust_factor <= 0) {
    stop("p_adjust_factor must be a single positive number.")
  }

  df <- df[
    is.finite(df$x_vec) & is.finite(df$y_vec) & !is.na(df$Groups),
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    stop("No complete x-y pairs with a valid group were found.")
  }

  correlation_label <- function(group_name) {
    group_df <- df[df$Groups == group_name, , drop = FALSE]

    if (nrow(group_df) < 3 ||
        is.na(sd(group_df$x_vec)) || sd(group_df$x_vec) == 0 ||
        is.na(sd(group_df$y_vec)) || sd(group_df$y_vec) == 0) {
      prefix <- paste0(group_name, ": r = NA, ")
      return(bquote(.(prefix) * italic(P) * " = NA"))
    }

    corr_result <- cor.test(
      group_df$x_vec,
      group_df$y_vec,
      method = "pearson"
    )
    r_value <- unname(corr_result$estimate)
    p_value <- corr_result$p.value
    sig_label <- if (p_value * p_adjust_factor < 0.05) {
      "**"
    } else if (p_value < 0.05) {
      "*"
    } else {
      ""
    }

    p_text <- if (p_value < 0.001) " < 0.001" else sprintf(" = %.3f", p_value)
    prefix <- sprintf("%s: r = %.3f%s, ", group_name, r_value, sig_label)
    bquote(.(prefix) * italic(P) * .(p_text))
  }

  mdd_label <- correlation_label("MDD")
  hc_label <- correlation_label("HC")
  correlation_labels <- as.expression(list(mdd_label, hc_label))
  stats_element <- element_text(
    size = stats_label_size, hjust = 0, margin = margin(b = 6)
  )
  class(stats_element) <- c("element_corr_stats", class(stats_element))

  p <- ggplot(
    data = df,
    mapping = aes(
      x = x_vec,
      y = y_vec,
      group = Groups,
      shape = Groups,
      color = Groups
    )
  ) +
    theme_classic(base_size = 30) +
    geom_point(size = 6) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(
      values = c("MDD" = "#E73F1E", "HC" = "#0D47A1")
    ) +
    ggtitle(title) +
    labs(
      x = xlab,
      y = ylab,
      subtitle = correlation_labels
    ) +
    theme(
      axis.title.x = element_text(size = x_axis_label_size),
      axis.title.y = element_text(size = y_axis_label_size),
      axis.text.x = element_text(size = x_tick_label_size),
      axis.text.y = element_text(size = y_tick_label_size),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = stats_element
    )
  
  if (!is.null(out_file)) {
    ggsave(out_file, plot = p, width = 25, height = 20, units = "cm", dpi = 300)
  }
  return(p)
}

plot_corr_txt_1group <- function(
    txt_file,
    out_file = NULL,
    title = "",
    xlab = "Scale score",
    ylab = "MRI metric",
    point_color = "#0D47A1",
    method = c("pearson", "spearman"),
    p_spin = NULL
) {

  library(ggplot2)

  df <- read.table(
    txt_file,
    header = FALSE,
    sep = "",
    col.names = c("x_vec", "y_vec")
  )

  if (ncol(df) != 2) {
    stop("The input text file must contain exactly two columns: x and y.")
  }

  df$x_vec <- as.numeric(df$x_vec)
  df$y_vec <- as.numeric(df$y_vec)
  df <- df[is.finite(df$x_vec) & is.finite(df$y_vec), , drop = FALSE]

  if (nrow(df) < 3) {
    stop("At least three complete x-y pairs are required.")
  }
  if (sd(df$x_vec) == 0 || sd(df$y_vec) == 0) {
    stop("Correlation cannot be calculated when x or y has zero variance.")
  }

  method <- match.arg(method)
  corr_result <- cor.test(df$x_vec, df$y_vec, method = method, exact = FALSE)
  r_value <- unname(corr_result$estimate)
  p_value <- corr_result$p.value

  r_label <- if (method == "spearman") {
    sprintf("\"Spearman's \"*rho == %.3f", r_value)
  } else {
    sprintf("italic(r) == %.3f", r_value)
  }
  # A supplied spatial-rotation p-value replaces the parametric test p-value.
  if (!is.null(p_spin)) {
    if (!is.numeric(p_spin) || length(p_spin) != 1L ||
        !is.finite(p_spin) || p_spin < 0 || p_spin > 1) {
      stop("p_spin must be one finite probability between 0 and 1.")
    }
    p_value <- p_spin
    p_label <- sprintf("italic(P)[spin] == %.3g", p_value)
  } else if (p_value < 0.001) {
    p_label <- "italic(p) < 0.001"
  } else {
    p_label <- sprintf("italic(p) == %.3f", p_value)
  }

  x_range <- range(df$x_vec)
  y_range <- range(df$y_vec)
  x_span <- diff(x_range)
  y_span <- diff(y_range)

  x_mid <- mean(x_range)
  y_mid <- mean(y_range)
  corner_counts <- c(
    top_left = sum(df$x_vec <= x_mid & df$y_vec >= y_mid),
    top_right = sum(df$x_vec > x_mid & df$y_vec >= y_mid),
    bottom_left = sum(df$x_vec <= x_mid & df$y_vec < y_mid),
    bottom_right = sum(df$x_vec > x_mid & df$y_vec < y_mid)
  )
  label_corner <- names(which.min(corner_counts))

  label_on_left <- grepl("left$", label_corner)
  label_on_top <- grepl("^top", label_corner)
  label_x <- if (label_on_left) {
    x_range[1] + 0.04 * x_span
  } else {
    x_range[2] - 0.04 * x_span
  }
  label_y <- if (label_on_top) {
    c(y_range[2] - 0.04 * y_span, y_range[2] - 0.13 * y_span)
  } else {
    c(y_range[1] + 0.13 * y_span, y_range[1] + 0.04 * y_span)
  }
  label_hjust <- if (label_on_left) 0 else 1

  p <- ggplot(
    data = df,
    mapping = aes(x = x_vec, y = y_vec)
  ) +
    theme_classic(base_size = 30) +
    geom_point(size = 6, color = point_color) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    annotate(
      "text",
      x = label_x,
      y = label_y[1],
      label = r_label,
      parse = TRUE,
      hjust = label_hjust,
      size = 8
    ) +
    annotate(
      "text",
      x = label_x,
      y = label_y[2],
      label = p_label,
      parse = TRUE,
      hjust = label_hjust,
      size = 8
    ) +
    ggtitle(title) +
    labs(
      x = xlab,
      y = ylab
    ) +
    theme(plot.title = element_text(hjust = 0.5))

  if (!is.null(out_file)) {
    ggsave(
      out_file,
      plot = p,
      width = 25,
      height = 20,
      units = "cm",
      dpi = 300
    )
  }

  return(p)
}

################### sEEG LMM t-statistic bar plot ###################
plot_seeg_lmm_tstat_bar <- function(
    df,
    out_file = NULL,
    bar_color = "#0D47A1",
    width = 120,
    height = 100,
    dpi = 600,
    title = NULL,
    y_limits = NULL,
    x_tick_label_size = 12,
    y_tick_label_size = 12,
    x_axis_label_size = 16,
    y_axis_label_size = 16,
    xlab = NULL,
    ylab = expression(italic(t) * " statistic")
) {
  
  library(ggplot2)
  library(dplyr)
  
  if (!is.data.frame(df) || nrow(df) == 0L) {
    stop("df must be a non-empty data frame.")
  }

  required_columns <- c("Band", "pValue", "pValue_FDR_BH")
  missing_columns <- setdiff(required_columns, colnames(df))
  
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing required columns: ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  if ("tStat" %in% colnames(df)) {
    t_column <- "tStat"
  } else if ("tStats" %in% colnames(df)) {
    t_column <- "tStats"
  } else {
    stop("Missing t-statistic column: expected 'tStat' or 'tStats'.")
  }
  
  band_order <- c(
    "delta",
    "theta",
    "alpha",
    "beta",
    "low gamma",
    "high gamma"
  )
  
  df <- df %>%
    mutate(
      Band = trimws(as.character(Band)),
      t_value = as.numeric(.data[[t_column]]),
      sig_label = case_when(
        pValue_FDR_BH < 0.05 ~ "**",
        pValue < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  if (any(!is.finite(df$t_value))) {
    stop("The t-statistic column contains missing or non-finite values.")
  }
  
  if (anyNA(df$Band) || any(!nzchar(df$Band)) || anyDuplicated(df$Band)) {
    stop("Band must contain one non-empty, unique value per row.")
  }
  if (!is.null(y_limits) &&
      (length(y_limits) != 2L || any(!is.finite(y_limits)) ||
       y_limits[1] >= y_limits[2] ||
       y_limits[1] > min(0, df$t_value) ||
       y_limits[2] < max(0, df$t_value))) {
    stop("y_limits must be two increasing finite values spanning zero and all bars.")
  }

  band_order <- c(
    intersect(band_order, df$Band),
    setdiff(unique(df$Band), band_order)
  )
  
  df <- df %>%
    mutate(
      Band = factor(Band, levels = band_order)
    )
  
  label_offset <- 0.06 * max(abs(c(df$t_value, y_limits)))
  if (!is.finite(label_offset) || label_offset == 0) {
    label_offset <- 0.1
  }
  
  p <- ggplot(
    df,
    aes(
      x = Band,
      y = t_value
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      colour = "black"
    ) +
    geom_col(
      width = 0.65,
      fill = bar_color,
      colour = NA
    ) +
    geom_text(
      aes(
        y = t_value + ifelse(
          t_value >= 0,
          label_offset,
          -label_offset
        ),
        label = sig_label
      ),
      vjust = ifelse(df$t_value >= 0, 0, 1),
      size = 6,
      fontface = "bold",
      family = "Arial"
    ) +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.12, 0.18))
    ) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    theme_classic(base_size = 9) +
    theme(
      text = element_text(
        family = "Arial",
        colour = "black"
      ),
      axis.line = element_line(
        linewidth = 0.5,
        colour = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.5,
        colour = "black"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = x_tick_label_size,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = y_tick_label_size,
        colour = "black"
      ),
      axis.title.x = element_text(
        size = x_axis_label_size,
        face = "bold",
        colour = "black"
      ),
      axis.title.y = element_text(
        size = y_axis_label_size,
        face = "bold",
        colour = "black"
      ),
      legend.position = "none",
      plot.margin = margin(8, 8, 5, 5)
    )
  
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(
      filename = out_file,
      plot = p,
      width = width,
      height = height,
      units = "mm",
      dpi = dpi,
      bg = "white"
    )
  }
  
  return(p)
}


################### ROI FC network t-statistic bar plot ###################
# Input: one site/group x cluster x GSR subset of the network_stats CSV,
# with network, t, p and p_fdr_7networks columns (Suzhou or Beijing).
# Significance follows plot_seeg_lmm_tstat_bar: ** FDR p < .05; * raw p < .05.
# Font sizes are in points; output width/height are in mm. Returns a ggplot.
# Example: roi_FC_network_bar(subset(stats, site == "IPCAS" & analysis ==
#   "positive_gt3p29_cluster1_n69_mask.func_FunSurfWCF"),
#   bar_color = c(Vis = "#781286", SomMot = "#4682B4", DorsAttn = "#00760E",
#                 SalVentAttn = "#C43AFA", Limbic = "#DCF8A4",
#                 Cont = "#E69422", Default = "#CD3E4E"))
roi_FC_network_bar <- function(
    df,
    out_file = NULL,
    bar_color = "#0D47A1",
    width = 120,
    height = 100,
    dpi = 600,
    title = NULL,
    y_limits = NULL,
    x_tick_label_size = 12,
    y_tick_label_size = 12,
    x_axis_label_size = 16,
    y_axis_label_size = 16,
    network_order = c("Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"),
    xlab = NULL,
    ylab = expression(italic(t) * " statistic")
) {
  
  library(ggplot2)
  library(dplyr)
  
  if (!is.data.frame(df) || nrow(df) == 0L) {
    stop("df must be a non-empty data frame.")
  }

  required_columns <- c("network", "t", "p", "p_fdr_7networks")
  missing_columns <- setdiff(required_columns, colnames(df))
  
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing required columns: ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  if (!is.numeric(df$t) || !is.numeric(df$p) ||
      !is.numeric(df$p_fdr_7networks)) {
    stop("t, p and p_fdr_7networks must be numeric.")
  }
  for (column in c("p", "p_fdr_7networks")) {
    values <- df[[column]]
    if (any(!is.na(values) & (!is.finite(values) | values < 0 | values > 1))) {
      stop(paste0(column, " must contain probabilities in [0, 1] or NA."))
    }
  }
  if (!is.character(network_order) || anyNA(network_order) ||
      any(!nzchar(network_order)) || anyDuplicated(network_order)) {
    stop("network_order must contain unique, non-empty network names.")
  }

  df <- df %>%
    mutate(
      network = trimws(as.character(network)),
      t_value = t,
      sig_label = case_when(
        p_fdr_7networks < 0.05 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  if (any(!is.finite(df$t_value))) {
    stop("The t-statistic column contains missing or non-finite values.")
  }
  
  if (anyNA(df$network) || any(!nzchar(df$network)) || anyDuplicated(df$network)) {
    stop("network must contain one non-empty, unique value per row; subset to one site/group and analysis first.")
  }
  if (!is.null(y_limits) &&
      (length(y_limits) != 2L || any(!is.finite(y_limits)) ||
       y_limits[1] >= y_limits[2] ||
       y_limits[1] > min(0, df$t_value) ||
       y_limits[2] < max(0, df$t_value))) {
    stop("y_limits must be two increasing finite values spanning zero and all bars.")
  }

  network_order <- c(
    intersect(network_order, df$network),
    setdiff(unique(df$network), network_order)
  )
  
  df <- df %>%
    mutate(
      network = factor(network, levels = network_order)
    )
  
  # A single colour fills all bars. Unnamed vectors follow network_order;
  # named vectors map colours by network, independently of input row order.
  plotted_networks <- levels(df$network)
  if (!is.character(bar_color) || !length(bar_color) || anyNA(bar_color)) {
    stop("bar_color must be a colour string or a vector of colours.")
  }
  if (!is.null(names(bar_color))) {
    if (anyNA(names(bar_color)) || any(!nzchar(names(bar_color))) ||
        anyDuplicated(names(bar_color)) ||
        !all(plotted_networks %in% names(bar_color))) {
      stop("Named bar_color must provide a unique colour for every plotted network.")
    }
    fill_colors <- bar_color[plotted_networks]
  } else if (length(bar_color) == 1L) {
    fill_colors <- setNames(rep(bar_color, length(plotted_networks)), plotted_networks)
  } else if (length(bar_color) == length(plotted_networks)) {
    fill_colors <- setNames(bar_color, plotted_networks)
  } else {
    stop("Unnamed bar_color must have length 1 or one colour per plotted network.")
  }
  tryCatch(grDevices::col2rgb(fill_colors), error = function(e) {
    stop("bar_color contains an invalid colour.", call. = FALSE)
  })

  label_offset <- 0.06 * max(abs(c(df$t_value, y_limits)))
  if (!is.finite(label_offset) || label_offset == 0) {
    label_offset <- 0.1
  }
  
  p <- ggplot(
    df,
    aes(
      x = network,
      y = t_value
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      colour = "black"
    ) +
    geom_col(
      width = 0.65,
      aes(fill = network),
      colour = NA
    ) +
    scale_fill_manual(values = fill_colors, guide = "none") +
    geom_text(
      aes(
        y = t_value + ifelse(
          t_value >= 0,
          label_offset,
          -label_offset
        ),
        label = sig_label
      ),
      vjust = ifelse(df$t_value >= 0, 0, 1),
      size = 6,
      fontface = "bold",
      family = "Arial"
    ) +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.12, 0.18))
    ) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    theme_classic(base_size = 9) +
    theme(
      text = element_text(
        family = "Arial",
        colour = "black"
      ),
      axis.line = element_line(
        linewidth = 0.5,
        colour = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.5,
        colour = "black"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = x_tick_label_size,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = y_tick_label_size,
        colour = "black"
      ),
      axis.title.x = element_text(
        size = x_axis_label_size,
        face = "bold",
        colour = "black"
      ),
      axis.title.y = element_text(
        size = y_axis_label_size,
        face = "bold",
        colour = "black"
      ),
      legend.position = "none",
      plot.margin = margin(8, 8, 5, 5)
    )
  
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(
      filename = out_file,
      plot = p,
      width = width,
      height = height,
      units = "mm",
      dpi = dpi,
      bg = "white"
    )
  }
  
  return(p)
}
