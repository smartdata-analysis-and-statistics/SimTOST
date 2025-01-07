#' @title Plot Power vs Sample Size for Simulation Results
#' @description Generates a detailed plot showing the relationship between power and total sample size for each comparator and the overall combined comparators.
#' The plot also includes confidence intervals for power estimates and highlights the target power with a dashed line for easy visual comparison.
#'
#' @param x An object of class `simss` containing simulation results.
#' @param \ldots Additional arguments to be passed to the `plot.simss` function for customization.
#'
#' @return A `ggplot` object illustrating:
#'   - Power (y-axis) vs. Total Sample Size (x-axis) for individual endpoints and comparators.
#'   - Error bars representing the 95% confidence interval of the power estimates.
#'   - A dashed horizontal line indicating the target power for comparison.
#'   - Faceted panels for each comparator, making it easy to compare results across different groups.
#'
#' @details
#' The plot dynamically adjusts to exclude unnecessary components, such as redundant endpoints or comparators with insufficient data, ensuring clarity and simplicity.
#' The `ggplot2` framework is used for visualizations, allowing further customization if needed.
#'
#' @author
#' Johanna Muñoz \email{johanna.munoz@fromdatatowisdom.com}
#'
#' @importFrom data.table .SD
#' @export plot.simss
#'
#' @export
plot.simss <- function(x, ...){
  qnam = n_iter = n_total = t_true = power = Endpoint = power_LCI = power_UCI = NULL # due to NSE notes in R CMD check

  table.test <- x$table.test
  nsim <- as.numeric(x$param.d["nsim"])
  tpower <- as.numeric(x$param.d["power"])
  # Calculate totaly test across all comparators= power
  qnam <- colnames(table.test)[grep("^totaly",colnames(table.test))]


  # Get a summary across all the n_iter with confidence interval power
  namexc <- c(colnames(table.test)[grep("^[^(mu_|sd_|eql_|equ_|n_)]",colnames(table.test))],"n_total")
  summary <- table.test[, lapply(.SD, FUN=function(x){sum(x, na.rm=TRUE)}), by= n_total ][,c(namexc),with=FALSE]
  plotdata <- data.table::melt(summary, id.vars = c("t_true","n_total"))
  powerfun <- function(x) {
    bin_test <- stats::prop.test(x = x, n = nsim, correct = TRUE)
    c(bin_test$estimate[[1]],bin_test$conf[1],bin_test$conf[2])
  }

  powerv <- do.call(rbind,lapply(plotdata$value, powerfun))
  colnames(powerv) <- c("power","power_LCI","power_UCI")
  plotdata <- cbind(plotdata,powerv)
  plotdata$Endpoint <- sub("(.*)Comp:.*", "\\1", plotdata$variable)
  plotdata$Comparator <- sub(".*Comp:", "", plotdata$variable)

  plotdata$Endpoint<-ifelse(plotdata$Endpoint=="totaly","Total",plotdata$Endpoint)
  plotdata$Comparator<-ifelse(plotdata$Comparator=="totaly","All comparators",plotdata$Comparator)


  table.uni <- data.frame(unclass(table(unique(plotdata[,c("Comparator","Endpoint")]))))
  table.uni$sum <- rowSums(table.uni)
  uniquey<-row.names(table.uni)[table.uni$sum<=2]
  uniquey<-uniquey[uniquey!="All comparators"]

  if(nrow(table.uni)==2){
    plotdata<-plotdata[ "Comparator" != "All comparators"]
  }

  if(length(uniquey)!=0){
    for(comp in uniquey){
      plotdata <- plotdata[!("Comparator"==comp & "Endpoint"=="Total")]
    }
  }


  plot <- ggplot2::ggplot(plotdata,ggplot2::aes(x=n_total,y=power,color=Endpoint))+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=power_LCI, ymax= power_UCI),width=0.1)+
    ggplot2::geom_point(size=0.3)+
    ggplot2::geom_line(linewidth=0.1)+
    ggplot2::geom_hline(yintercept=tpower,linetype="dashed")+
    ggplot2::facet_grid(.~Comparator)+
    ggplot2::theme_minimal()+
    ggplot2::xlab("Total sample size")+ ggplot2::ylab("Power")

  if(length(unique(plotdata$Endpoint))>1){
    plot <- plot+ggplot2::theme(legend.position="bottom")
  }else{
    plot <- plot+ggplot2::theme(legend.position = "none")
  }

  return(plot)
}
