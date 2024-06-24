#' Setting up the package dartR.sexlinked
#'
#' Setting up dartR.sexlinked
#' @keywords internal
#' @importFrom utils packageVersion
#' @importFrom stats complete.cases
#' @importFrom methods getPackageName
#' @import ggplot2
#' @import adegenet
#' @import dartR.base
#' @import dartR.data



#needed to avoid error
zzz <- NULL

# defining function %dopar% from foreach package to pass CRAN checks
`%dopar%` <- foreach::`%dopar%`

build = "Jody"
error <- crayon::red
warn <- crayon::yellow
report <- crayon::green
important <- crayon::blue
code <- crayon::cyan


# WELCOME MESSAGE
.onAttach <- function(...) {
  pn <- getPackageName()
  packageStartupMessage(important(
    paste(
      "**** Welcome to",pn,"[Version",
      packageVersion(pn),
      "] ****\n"
    )
  ))
}

