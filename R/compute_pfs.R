#' Compute progression-free-survival rate given sample
#'
#' `compute_pfs()` computes the progression-free-survival rate at specified
#' times given a paramter sample.
#'
#' @template param-model
#' @param t a vector of time-points at which the PFS rate should be computed
#' @template param-parameter_sample
#' @template param-warmup
#' @template param-nsim
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return a data frame with samples of PFS rates at each of the time points
#' in the vector t.
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' smpl <- sample_prior(mdl, nsim = 500, seed = 34L)
#' dplyr::filter(
#'   compute_pfs(mdl, t = seq(0, 12), parameter_sample = smpl),
#'   iter == 1
#' )
#'
#' @export
compute_pfs <- function(
  model,
  t,
  parameter_sample = NULL,
  warmup = 500L,
  nsim = 1000L,
  seed = NULL,
  ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  if (is.null(parameter_sample)) {
    parameter_sample <- sample_prior(model,
                                     warmup = warmup, nsim = nsim, seed = seed,
                                     ...)
  }
  group_ids <- model$group_id
  n_groups <- length(group_ids)
  n_smpl <- dim(parameter_sample)[[1]]
  n_t <- length(t)
  res <- array(0.0, dim = c(n_smpl, n_t, n_groups))
  p <- rstan::extract(parameter_sample, "p")[[1]]
  shape <- rstan::extract(parameter_sample, "shape")[[1]]
  scale <- rstan::extract(parameter_sample, "scale")[[1]]
  for (i in 1:n_smpl) {
    for (j in 1:n_groups) {
      # calculate pfs for i-th sample, j-th group and add to result
      res[i, , j] <- pfs(t, p[i, j], shape[i, j, ], scale[i, j, ])
    }
  }
  # unpack 3d array into tibble
  tbl_pfs <- tibble(
    iter = integer(),
    group_id = character(),
    t = numeric(),
    pfs = numeric()
  )
  dimnames(res) <- list(1:n_smpl, 1:n_t, group_ids)
  for (i in seq_along(group_ids)) {
    tbl_pfs <- dplyr::bind_rows(
      tbl_pfs,
      as_tibble(t(res[, , i])) %>%
        mutate(t = t) %>%
        tidyr::pivot_longer(
          -t, names_to = "iter", values_to = "pfs"
        ) %>%
        mutate(group_id = group_ids[i], iter = as.integer(.data$iter))
    )
  }
  return(tbl_pfs)
}
