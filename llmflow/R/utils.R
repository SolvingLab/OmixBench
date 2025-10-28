.onAttach <- function(libname, pkgname) {
  suppressPackageStartupMessages(library(ellmer))
  msg <- paste0(
    "===================================================================\n",
    "========================= Welcome to llmflow ======================\n",
    "===================================================================\n"
  )
  base::packageStartupMessage(msg)
}
