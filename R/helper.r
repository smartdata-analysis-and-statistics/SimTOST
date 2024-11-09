#' Display the summary results of the sample size estimation
#' @details
#' This function displays the summary results of the sample size estimation.
#'
#'
#' @param x An object of class "simss"
#' @param ...  Optional arguments to be passed from or to other methods.
#' @author
#' Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @method print simss
#' @export
print.simss <- function(x, ...) {

  alpha <- x$param.d[["alpha"]]
  tpower <- x$param.d[["power"]]
  power <- round(x$response[["power"]], digits = 3)
  lpower <- round(x$response[["power_LCI"]], digits = 3)
  upower <- round(x$response[["power_UCI"]], digits = 3)
  sst <- x$response[["n_total"]]

  message(cat(paste0("Given a ",100*tpower,"%  target power with 100(1-2*",alpha, ")% confidence level.")))
  message(cat(paste0("The total required sample size to achieve ",100*power,"% power is ",sst," sample units.\n")))
  print(x$response[,-1], row.names = FALSE)
}


