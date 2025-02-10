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

  cat("Sample Size Calculation Results\n")
  cat("-------------------------------------------------------------\n")
  cat(paste0("Study Design: ", x$param.d$dtype, " trial targeting ",100*tpower,"% power with a ",100*alpha, "% type-I error.\n\n"))
  cat("Comparisons:\n")
  for (i in 1:length(x$param$list_comparator)) {
    nendp <- length(x$param$list_y_comparator[[i]])
    cat("  ", paste(x$param$list_comparator[[i]], collapse = " vs. "), "\n")
    cat("    - Endpoints Tested:", paste(x$param$list_y_comparator[[i]], collapse = ", "), "\n")
    if (nendp > 1 & x$param.d$k[i] < nendp) {
      cat("      (multiple primary endpoints, k = ", x$param.d$k[i], ")\n")
    } else if (nendp > 1 & x$param.d$k[i] == nendp) {
      if (x$param.d$adjust == "seq") {
        cat("      (hierarchical testing applied to multiple endpoints, m =", x$param.d$k[i], ")\n")
      } else {
        cat("      (multiple co-primary endpoints, m = ", x$param.d$k[i], ")\n")
      }
    }
    if (x$param.d$adjust == "no") {
      multiplicity_correction <- "No Adjustment"
      alphau <- x$param.d$alpha
    } else if (x$param.d$adjust == "bon") {
      multiplicity_correction <- "Bonferroni"
      alphau <- x$param.d$alpha/nendp
    } else if (x$param.d$adjust == "sid") {
      multiplicity_correction <- "Sidak"
      alphau <- 1 - (1 - x$param.d$alpha)^(1/nendp)
    } else if (x$param.d$adjust == "k") {
      multiplicity_correction <- "k-adjustment"
      alphau <- x$param.d$k[i]*x$param.d$alpha/nendp
    } else if (x$param.d$adjust == "seq") {
      multiplicity_correction <- "Sequential"
      alphau <- x$param.d$alpha*x$param$weight_seq
    }

    # Generate the description of the multiplicity correction in case of multiple primary endpoints
    if (nendp > 1 & (x$param.d$k[i] < nendp | x$param.d$adjust == "seq" )) {
      cat("    - Multiplicity Correction:", multiplicity_correction, "\n")
      cat("      Adjusted Significance Levels: alpha =", paste(format(alphau, digits = 3, nsmall = 3), collapse = "; "), "\n\n")
    }

  }



  cat("-------------------------------------------------------------\n")
  print(output, row.names = FALSE)  # Suppress row numbers
  cat("-------------------------------------------------------------\n")
  #cat("-------------------------------------------------------------\n")
  #print(x$response[,-1], row.names = FALSE)
}
