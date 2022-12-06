#' Print an srpmodel
#'
#' @param x model to print
#' @template param-dotdotdot
#' @examples
#' print(create_srpmodel(A = define_srp_prior()))
#' @export
print.srpmodel <- function(x, ...) cat(format(x, ...), "\n") # nocov

#' @examples
#' format(create_srpmodel(A = define_srp_prior()))
#' @rdname print.srpmodel
#' @export
format.srpmodel <- function(x, ...) {
  sprintf("srpmodel<%s>", paste(x$group_id, collapse = ","))
}
