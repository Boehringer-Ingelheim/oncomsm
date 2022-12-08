#' Print an srpmodel
#'
#' @param x model to print
#' @template param-dotdotdot
#'
#' @return `format()` returns a character string representation of the object,
#' `print()` prints to the console and returns the object itself invisibly.
#'
#' @examples
#' print(create_srpmodel(A = define_srp_prior()))
#' @export
print.srpmodel <- function(x, ...) {
  cat(format(x, ...), "\n") # nocov
  return(invisible(x)) # nocov
}

#' @examples
#' format(create_srpmodel(A = define_srp_prior()))
#' @rdname print.srpmodel
#' @export
format.srpmodel <- function(x, ...) {
  sprintf("srpmodel<%s>", paste(x$group_id, collapse = ","))
}
