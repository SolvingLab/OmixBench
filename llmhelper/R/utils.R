#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

#' @export
`%||%` <- function(x, y) if (is.null(x)) y else x
