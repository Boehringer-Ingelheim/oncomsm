#' @export
plot_visits <- function(tbl_data, epsilon = 5/30) {
  t_max <- max(tbl_data$t)
  res <- tbl_data %>%
    group_by(subject_id) %>%
    mutate(
      t_end = case_when(
        row_number() == n() & !eof ~ t_max + epsilon,
        row_number() == n() ~ t + epsilon,
        row_number() != n() ~ lead(t)
      )
    ) %>%
    ungroup() %>%
    rename(t_start = t) %>%
    select(subject_id, group_id, t_start, t_end, status, eof)
  res <- res %>%
    group_by(group_id) %>%
    mutate(
      subject_id = as.numeric(factor(subject_id, levels = unique(subject_id[order(t_start)]))),
      group_id = factor(group_id)
    )
  ggplot(res) +
    geom_rect(
      aes(
        xmin = t_start, xmax = t_end,
        ymin = as.numeric(subject_id) - .5,
        ymax = as.numeric(subject_id) + .5,
        fill = status
      )
    ) +
    scale_fill_manual(
      "Status",
      values = c(R = "darkgreen", S = "gray", P = "red"),
      na.translate = FALSE
    ) +
    labs(
      x = "months since start of trial",
      y = "subject"
    ) +
    facet_wrap(~group_id, scales = "free_y") +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.spacing = unit(.2, "lines"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}
