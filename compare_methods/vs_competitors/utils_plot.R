
time_plots = function(df, baseline = "BRIDGES") {
  p_time = res %>%
    dplyr::filter(metric == "seconds") %>%
    ggplot(mapping = aes(x = as.factor(ncells), y = value, col = algorithm)) +
    #geom_jitter() +
    geom_beeswarm() +
    # geom_quasirandom() +
    #geom_line() +
    scale_y_continuous(transform = "log10") +
    theme_bw() +
    scale_colour_nejm() +
    labs(x = "Number of cells", y = "Time (s)", col = "Algorithm")
  p_time

  p_time_frac = res %>%
    dplyr::filter(metric == "seconds") %>%
    dplyr::group_by(sim_id) %>%
    dplyr::mutate(f = value / value[algorithm == baseline]) %>%
    ggplot(mapping = aes(x = as.factor(ncells), y = f, col = algorithm)) +
    #geom_jitter() +
    geom_beeswarm() +
    # geom_quasirandom() +
    #geom_line() +
    scale_y_continuous(transform = "log10") +
    theme_bw() +
    scale_colour_nejm() +
    labs(x = "Number of cells", y = paste0("Time-fold-change w.r.t ", baseline), col = "Algorithm")

  list(p_time = p_time, p_time_frac = p_time_frac)
}


metric_plot = function(df, metric_name) {
  df$bfb_prop = df$bfb_rate / (df$normal_rate + df$del_rate + df$amp_rate + df$bfb_rate)
  df$BFB = ifelse(df$bfb_prop > .01, "High", "Low")

  p_norm = df %>%
    dplyr::filter(metric == metric_name) %>%
    dplyr::group_by(sim_id) %>%
    dplyr::mutate(normed_metric = value / max(value)) %>%
    ggplot(mapping = aes(x = algorithm, y = normed_metric, col = algorithm)) +
    geom_boxplot() +
    geom_jitter(aes(col=algorithm), alpha = 1, size = .8) +
    theme_bw() +
    labs(x = "", y=paste0("Norm'd ", metric_name), col = "Algorithm") +
    theme(legend.position = "bottom") +
    facet_wrap(~paste0("BFB ", BFB)) +
    scale_colour_nejm() +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  p_norm

  p_absolute = df %>%
    dplyr::filter(metric == metric_name) %>%
    ggplot(mapping = aes(x = factor(ncells), y = value, col = algorithm)) +
    geom_boxplot() +
    ggh4x::facet_nested(~"BFB proportion"+bfb_prop) +
    #facet_wrap(~bfb_prop) +
    scale_colour_nejm() +
    theme_bw() +
    labs(x = "Number of cells", y=metric_name, col = "Algorithm") +
    theme(legend.position = "bottom") +
    scale_colour_nejm()
  p_absolute

  list(p_norm = p_norm, p_absolute = p_absolute)
}


plot_trend = function(df, metric_name, a = .25, logy = F, logx = F)  {
  df$bfb_prop = df$bfb_rate / (df$normal_rate + df$del_rate + df$amp_rate + df$bfb_rate)
  p = df %>%
    dplyr::filter(metric == metric_name) %>%
    dplyr::group_by(algorithm, ncells, bfb_prop, metric, normal_rate) %>%
    na.omit() %>%
    dplyr::summarise(y = mean(value), yl = quantile(value, a, na.rm = TRUE), yh = quantile(value, 1-a, na.rm = TRUE)) %>%
    ggplot(mapping = aes(x = ncells, y=y, ymax=yh, ymin=yl, col=algorithm)) +
    geom_pointrange() +
    geom_line() +
    #scale_x_continuous(transform = "log10") +
    ggh4x::facet_nested(~"BFB prop"+bfb_prop, scales = "free_y") +
    theme_bw() +
    labs(x = "Number of leaves", y=metric_name, col = "Method") +
    theme(legend.position = "bottom")

  if (logy) p = p + scale_y_continuous(transform = "log10")
  if (logx) p = p + scale_x_continuous(transform = "log10")
  p
}

plot_normed_boxplots = function(df, metric_name, split_by = NULL) {
  df$bfb_prop = df$bfb_rate / (df$normal_rate + df$del_rate + df$amp_rate + df$bfb_rate)
  df$BFB = ifelse(df$bfb_prop > .01, "High", "Low")

  p = df %>%
    dplyr::filter(metric == metric_name) %>%
    dplyr::group_by(sim_id) %>%
    dplyr::mutate(normed_metric = value / max(value)) %>%
    ggplot(mapping = aes(x = algorithm, y = normed_metric)) +
    geom_boxplot(col = "black") +
    geom_jitter(aes(col=algorithm), alpha = 1, size = .8) +
    theme_bw() +
    labs(x = "", y=paste0("Norm'd ", metric_name), col = "Method") +
    theme(legend.position = "bottom")

  if (!is.null(split_by)) p = p + facet_wrap(~paste(split_by, .data[[split_by]]))
  p
}
