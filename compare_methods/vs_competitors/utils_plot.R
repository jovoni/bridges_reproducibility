

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
