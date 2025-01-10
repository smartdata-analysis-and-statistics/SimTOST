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
  cat(paste0("Study Design: ", x$param.d$dtype, " trial targeting ",100*tpower,"% power with a ",100*alpha, "% type-I error.\n"))
  cat("Endpoints Tested:\n")
  for (i in 1:length(x$param$list_comparator)) {
    nendp <- length(x$param$list_y_comparator[[i]])
    str <- paste0("  - ", paste(x$param$list_comparator[[i]], collapse = " vs. "), ": ",
                  paste(x$param$list_y_comparator[[i]], collapse = ", "))
    if (nendp > 1 & x$param.d$k[i] < nendp) {
      str <- paste0(str, " (multiple primary endpoints, k = ", x$param.d$k[i], ")\n")
    } else if (nendp > 1 & x$param.d$k[i] == nendp) {
      str <- paste0(str, " (multiple co-primary endpoints, m = ", x$param.d$k[i], ")\n")
    }
  }

  cat("Multiplicity Correction:\n")
  for (i in 1:length(x$param$list_comparator)) {
    nendp <- length(x$param$list_y_comparator[[i]])
    str <- paste0("  - ", paste(x$param$list_comparator[[i]], collapse = " vs. "), ": ")
    if (x$param.d$adjust == "no") {
      str <- paste0(str, " no adjustment (alpha = ", x$param.d$alpha, ")\n")
    } else if (x$param.d$adjust == "bon") {
      alphau <- x$param.d$alpha/nendp
      str <- paste0(str, " bonferroni (alpha = ", format(alphau, digits = 3, nsmall = 3), ")\n")
    } else if (x$param.d$adjust == "sid") {
      alphau <- 1-(1-x$param.d$alpha)^(1/nendp)
      str <- paste0(str, " sidak (alpha = ", format(alphau, digits = 3, nsmall = 3), ")\n")
    } else if (x$param.d$adjust == "k") {
      alphau <- x$param.d$k[i]*x$param.d$alpha/nendp
      str <- paste0(str, " k-adjustment (alpha = ", format(alphau, digits = 3, nsmall = 3), ")\n")
    } else if (x$param.d$adjust == "seq") {
      str <- paste0(str, " sequential (alpha = ",  paste(x$param.d$alpha*x$param$weight_seq, collapse = "; "), ")\n")
    }
  }

  cat("-------------------------------------------------------------\n")
  print(output, row.names = FALSE)  # Suppress row numbers
  #cat("-------------------------------------------------------------\n")
  #print(x$response[,-1], row.names = FALSE)
}


generate_correction_description <- function(adjust, k = NULL, m = NULL, alpha = 0.05) {
  # Validate input
  if (!adjust %in% c("k", "bon", "sid", "no", "seq")) {
    stop("Invalid 'adjust' method. Must be one of: 'k', 'bon', 'sid', 'no', 'seq'.")
  }

  # Generate descriptions based on the adjustment method
  description <- switch(adjust,
                        "k" = {
                          if (is.null(k) || is.null(m)) stop("'k' and 'm' must be provided for 'k' adjustment.")
                          paste0("k-adjustment: Adjusted significance level based on the number of endpoints required for equivalence ",
                                 "(k = ", k, ") and the total number of tests (m = ", m, ") using the formula α_k = (k * α) / m. ",
                                 "Adjusted α = ", round((k * alpha) / m, 4), ".")
                        },
                        "bon" = {
                          if (is.null(m)) stop("'m' must be provided for Bonferroni adjustment.")
                          paste0("Bonferroni adjustment: Significance level adjusted to control the family-wise error rate (FWER) using ",
                                 "α/m, where m = ", m, ". Adjusted α = ", round(alpha / m, 4), ".")
                        },
                        "sid" = {
                          if (is.null(m)) stop("'m' must be provided for Sidak adjustment.")
                          sid_alpha <- 1 - (1 - alpha)^(1 / m)
                          paste0("Sidak adjustment: Adjusted significance level to control FWER under the assumption of independent tests. ",
                                 "Adjusted α = 1 - (1 - α)^(1/m), where m = ", m, ". Adjusted α = ", round(sid_alpha, 4), ".")
                        },
                        "no" = "No adjustment: No correction applied to the significance level.",
                        "seq" = {
                          paste0("Sequential adjustment: Tests are performed in a predefined sequence, starting with primary endpoints ",
                                 "and proceeding to secondary endpoints only if the primary endpoints meet the significance criteria. ",
                                 "Separate α adjustments are applied to each group of endpoints.")
                        }
  )

  return(description)
}



