# library(tidyverse)
#
# tbl_recist_data <- haven::read_sas("adrs_spd.sas7bdat") %>%
#   group_by(TRTA) %>%
#   filter(
#     PARAM == "Overall Response",
#     str_detect(TRTA, "COHORT"),
#     AVALC != "NE", # non evaluable
#     n() >= 20 # filter very small cohorts
#   ) %>% #pull(TRTA) %>% unique %>%  str_extract("(?<=COHORT )[A-Z]{1}")
#   ungroup() %>%
#   transmute(
#     subject_id = SUBJID,
#     group_id = str_extract(TRTA, "(?<=COHORT )[A-Z]{1}"),
#     visit_t = difftime(ADT, min(TRTSDT), units = "days") %>% {as.numeric(.)/30},
#     status = AVALC,
#     end_of_treatment = ONGFL != "Y"
#   )
#
# tbl_bla <- tbl_recist_data %>%
#   group_by(subject_id) %>%
#   summarize(tmin = min(visit_t), tmax = max(visit_t))
#
# tbl_recist_data %>%
#   filter(group_id == "A") %>%
#   ggplot() +
#     geom_segment(
#       aes(x = tmin, xend = tmax, y = subject_id, yend = subject_id), data = tbl_bla
#     ) +
#     geom_point(
#       aes(x = visit_t, y = subject_id, color = status, shape = end_of_treatment),
#       size = 5
#     ) +
#     labs(
#       x = "months since start of trial",
#       y = "subject",
#       caption = "Diamond indicates end-of-treatment."
#     ) +
#     scale_color_manual(
#       "RECIST",
#       values = c(
#         CR = "darkgreen",
#         PR = "blue",
#         SD = "orange",
#         PD = "red"
#       ),
#       na.translate = FALSE
#     ) +
#     facet_wrap(~group_id, ncol = 1) +
#     theme_bw() +
#     theme(
#       legend.position = "top"
#     )
#
# ggsave("tmp.pdf", width=8, height = 20)
