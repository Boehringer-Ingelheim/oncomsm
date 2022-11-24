library(tidyverse)
library(oncomsm)

group_ids <- c("B", "A")

tbl_data <- tibble(
  subject_id = as.character(c(1, 1, 1, 2, 2, 3, 3, 3, 4, 5, 5)),
  group_id = factor(c(rep("B", 5), rep("A", 6)), levels = group_ids),
  t = c(0, 1.2, 2.4, 0.3, 1.5, 0.7, 1.9, 3.1, 0, 1.2, 1.7),
  state = c("stable", "response", "EOF", "stable", "response", "stable", "stable", "progression", "stable", "stable", "response")
)

p <- c(0.1, 0.9)
shape <- matrix(c(
    1, 1, 1,
    1, 1, 1
  ), nrow = 2, byrow = TRUE)
scale <- matrix(c(
    3, 3, 5,
    3, 3, 5
  ), nrow = 2, byrow = TRUE)
visit_spacing <- c(1.2, 1.2)

# tbl_data

tmp <- tbl_data %>%
  oncomsm:::f(
    p, shape, scale, visit_spacing, 12*10
  ) %>%
  as_tibble() %>%
  mutate(
    group_id = as.character(group_id)
  )

mdl <- create_srp_model(
  c("B", "A"),
  logodds_mean = c(0, 0),
  median_time_to_next_event_mean = matrix(c(
    3, 3, 5,
    3, 3, 5
  ), nrow = 2, byrow = TRUE),
  visit_spacing = c(1.2, 1.2)
)

visits_to_mstate(tmp, mdl)
