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

  output <- data.frame(Parameter = c("Total Sample Size", "Achieved Power", "Power Confidence Interval"),
                       Value = c(sst, 100*power, paste0(100*lpower, " - ", 100*upower)))

  message("Sample Size Calculation Results")
  cat("-------------------------------------------------------------\n")
  cat(paste0("Study Design: ", x$param.d$dtype, " trial targeting ",100*tpower,"% power with a ",100*alpha, "% type-I error.\n"))
  cat("-------------------------------------------------------------\n")
  print(output, row.names = FALSE)  # Suppress row numbers
  #cat("-------------------------------------------------------------\n")
  #print(x$response[,-1], row.names = FALSE)
}



