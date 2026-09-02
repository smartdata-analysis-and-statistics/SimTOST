test_that("stability and Monte Carlo plots facet a sample-size curve", {
  skip_if_not_installed("ggplot2")
  curve <- simPower(
    n = c(10, 12), distribution = "norm", ctype = "DOM",
    mu_list = list(T = c(y1 = 1.1), R = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.5)),
    list_uequi.tol = list(T_vs_R = c(y1 = .5)),
    nsim = 10, ncores = 1, seed = 1
  )
  stability <- plot_stability(curve)
  mc_error <- plot_mc_error(curve)
  power_plot <- plot(curve)
  expect_gt(nrow(ggplot2::ggplot_build(stability)$layout$layout), 1)
  expect_gt(nrow(ggplot2::ggplot_build(mc_error)$layout$layout), 1)
  expect_s3_class(mc_error, "ggplot")
  expect_identical(unique(as.character(stability$data$Endpoint)), "Total")
  expect_identical(unique(as.character(mc_error$data$Endpoint)), "Total")
  expect_true(all(c("y1", "Total") %in%
                  unique(as.character(plot_stability(curve, endpoint = "all")$data$Endpoint))))
  expect_error(plot_stability(curve, endpoint = "missing"), "Unknown endpoint")
  expect_null(power_plot$labels$title)
  expect_null(power_plot$labels$subtitle)
  expect_match(mc_error$labels$title, "Achieved Monte Carlo error")
  expect_true(any(grepl("%", mc_error$data$label)))
  expect_match(stability$labels$subtitle, "Compared total sample sizes")
  expect_match(mc_error$labels$subtitle, "Compared total sample sizes")
  expect_error(plot_mc_error(curve, target_error = 0.025), "no longer available")
  expect_error(plot_stability(curve, overall = TRUE), "requires an 'All comparators'")
  expect_error(plot_mc_error(curve, overall = TRUE), "requires an 'All comparators'")
  expect_error(plot_stability(curve, overall = 1), "'overall'")
})

test_that("diagnostic plots compare objects with different nsim values", {
  make_power <- function(nsim) suppressWarnings(simPower(
    n = 10, distribution = "norm", ctype = "DOM",
    mu_list = list(T = c(y1 = 1.1), R = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.5)),
    list_uequi.tol = list(T_vs_R = c(y1 = .5)),
    nsim = nsim, ncores = 1, seed = nsim
  ))
  stability <- plot_stability(list(make_power(8), make_power(12)))
  mc_error <- plot_mc_error(list(make_power(8), make_power(12)))
  expect_gt(nrow(ggplot2::ggplot_build(stability)$layout$layout), 1)
  expect_gt(nrow(ggplot2::ggplot_build(mc_error)$layout$layout), 1)
  expect_error(plot_stability(list(make_power(8), "not a simulation")),
               "must contain")
})

test_that("simss diagnostics retain only the suggested sample size", {
  object <- structure(list(
    response = data.frame(n_total = 20),
    table.iter = data.frame(n_total = c(10, 20)),
    table.test = data.frame(
      n_total = c(10, 20, 20, 20),
      `totalyComp:T vs R` = c(0, 0, 1, 1),
      totaly = c(0, 0, 1, 1)
    ),
    param.d = list(nsim = 4)
  ), class = "simss")
  stability <- plot_stability(object)
  mc_error <- plot_mc_error(object)
  expect_identical(unique(as.character(stability$data$Endpoint)), "Total")
  expect_identical(unique(as.character(mc_error$data$Endpoint)), "Total")
  expect_match(stability$labels$subtitle, "Suggested total sample size = 20")
  expect_match(mc_error$labels$subtitle, "Suggested total sample size = 20")
  expect_equal(nrow(ggplot2::ggplot_build(stability)$layout$layout), 1)
})

test_that("continuous sample-size plots show intervals on every point", {
  skip_if_not_installed("ggplot2")
  table_test <- data.frame(
    t_true = rep(0, 4),
    n_total = c(10, 10, 20, 20),
    totaly = c(0, 1, 1, 1),
    check.names = FALSE
  )
  table_test[["totalyComp:T vs R"]] <- c(0, 1, 1, 1)
  table_test[["y1Comp:T vs R"]] <- c(0, 1, 1, 1)
  object <- structure(list(
    response = data.frame(n_total = 20),
    table.test = data.table::as.data.table(table_test),
    param.d = list(nsim = 2),
    param = list(list_comparator = list(c("T", "R")))
  ), class = "simss")

  built <- ggplot2::ggplot_build(plot(object, all = FALSE))
  expect_equal(nrow(built$layout$layout), 1L)
  overall_built <- ggplot2::ggplot_build(plot(object))
  expect_false("Comparator" %in% names(overall_built$layout$layout))
  expect_identical(unique(as.character(overall_built$plot$data$Endpoint)), "Total")
  point_layer <- built$data[[3L]]
  selected_layer <- built$data[[4L]]
  interval_layer <- built$data[[5L]]
  expect_true(all(c("ymin", "ymax") %in% names(interval_layer)))
  expect_equal(nrow(interval_layer), nrow(point_layer))
  expect_true(all(interval_layer$linewidth >= 2))
  expect_true(all(interval_layer$colour == "black"))
  expect_identical(built$plot$theme$legend.position, "none")
  expect_true(all(interval_layer$width >= 0.4))
  expect_equal(nrow(selected_layer), 2L)
  expect_null(built$plot$labels$title)
  expect_null(built$plot$labels$subtitle)
  expect_error(plot(object, target_power = 1.1), "target_power")
  expect_error(plot(object, display = "missing"), "Unknown comparator")
  expect_error(plot(object, endpoint = "missing"), "Unknown endpoint")
})

test_that("continuous sample-size plots sort candidates and highlight the first target crossing", {
  skip_if_not_installed("ggplot2")
  table_test <- data.frame(
    n_iter = c(20, 20, 10, 10, 30, 30),
    n_total = c(20, 20, 10, 10, 30, 30),
    t_true = c(2, 2, 0, 0, 2, 2),
    `totalyComp:T vs R` = c(1, 1, 0, 0, 1, 1),
    `y1Comp:T vs R` = c(1, 1, 0, 0, 1, 1),
    check.names = FALSE
  )
  object <- structure(list(
    # Deliberately point response at 30: the plot should derive the first
    # evaluated candidate reaching the target from the global total curve.
    response = data.frame(n_total = 30),
    table.test = data.table::as.data.table(table_test),
    param.d = list(nsim = 2),
    param = list(list_comparator = list(c("T", "R")))
  ), class = "simss")

  built <- ggplot2::ggplot_build(plot(object, target_power = .8))
  line_layer <- built$data[[2L]]
  selected_layer <- built$data[[4L]]
  expect_true(all(vapply(split(line_layer$x, line_layer$group),
                         function(values) all(diff(values) >= 0),
                         logical(1))))
  expect_equal(unique(selected_layer$x), 20)
})

test_that("sample-size plots use the stored planning target when omitted", {
  skip_if_not_installed("ggplot2")
  table_test <- data.frame(
    n_iter = rep(c(20, 30), each = 5),
    n_total = rep(c(20, 30), each = 5),
    t_true = c(rep(c(1, 1, 1, 1, 0), 1), rep(1, 5)),
    `totalyComp:T vs R` = c(rep(c(1, 1, 1, 1, 0), 1), rep(1, 5)),
    `y1Comp:T vs R` = c(rep(c(1, 1, 1, 1, 0), 1), rep(1, 5)),
    check.names = FALSE
  )
  object <- structure(list(
    response = data.frame(n_total = 30),
    table.test = data.table::as.data.table(table_test),
    param.d = list(nsim = 5, power = .9),
    param = list(list_comparator = list(c("T", "R")))
  ), class = "simss")

  inferred <- ggplot2::ggplot_build(plot(object))
  explicit <- ggplot2::ggplot_build(plot(object, target_power = .8))
  expect_equal(inferred$data[[1]]$yintercept, .9)
  expect_equal(unique(inferred$data[[4]]$x), 30)
  expect_equal(explicit$data[[1]]$yintercept, .8)
  expect_equal(unique(explicit$data[[4]]$x), 20)
})

test_that("power plot selection falls back when no plotted candidate reaches target", {
  data <- data.frame(
    n = c(20, 10), power = c(.6, .7), Endpoint = c("Total", "Total"),
    Comparator = c("All comparators", "All comparators")
  )
  expect_equal(
    SimTOST:::.select_power_plot_n(data, target_power = .8,
                                   n_col = "n", fallback = 20),
    20
  )
  expect_equal(
    SimTOST:::.select_power_plot_n(data[0, , drop = FALSE], target_power = .8,
                                   n_col = "n", fallback = 20),
    20
  )
})

test_that("sample-size selection uses the global count-comparison series", {
  data <- data.frame(
    n = c(20, 30, 20, 30),
    power = c(.95, .99, .70, .85),
    Endpoint = rep("Total", 4),
    Comparator = c("EMA", "EMA", "All comparisons", "All comparisons")
  )
  expect_equal(
    SimTOST:::.select_power_plot_n(
      data, target_power = .8, n_col = "n", fallback = 20
    ),
    30
  )
})

test_that("count sample-size plots highlight the global target crossing", {
  skip_if_not_installed("ggplot2")
  table_test <- data.frame(
    n_total = rep(c(20, 30), each = 4),
    `totalyComp:EMA` = 1,
    `totalyComp:US` = 1,
    `y1Comp:EMA` = 1,
    `y1Comp:US` = 1,
    totaly = c(1, 1, 0, 0, rep(1, 4)),
    check.names = FALSE
  )
  object <- structure(
    list(
      n_total = 30, target_power = .8,
      table.test = table_test
    ),
    class = "countss"
  )
  built <- ggplot2::ggplot_build(plot(object))
  expect_true(all(c("ymin", "ymax") %in% names(built$data[[4L]])))
  expect_true(all(built$data[[4L]]$linewidth >= 0.8))
  expect_equal(unique(built$data[[5]]$x), 30)
  expect_null(built$plot$labels$title)
  expect_null(built$plot$labels$subtitle)
})

test_that("distribution plots can display reconstructed TOST t-values", {
  object <- simPower(
    n = 12, distribution = "norm", ctype = "DOM",
    mu_list = list(T = c(y1 = 1.1), R = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.5)),
    list_uequi.tol = list(T_vs_R = c(y1 = .5)),
    nsim = 20, ncores = 1, seed = 1, keep_sim_data = TRUE
  )
  t_plot <- plot_distribution(object, estimand = "t_value",
                              n_trials = 10,
                              seed = 2026)
  expect_s3_class(t_plot, "ggplot")
  expect_setequal(unique(as.character(t_plot$data$boundary)), c("lower", "upper"))
  expect_match(t_plot$labels$x, "TOST t-statistic")
  expect_s3_class(plot_distribution(object, estimand = "t_value", arms = "T"),
                  "ggplot")
  expect_error(plot_distribution(object, overall = TRUE),
               "not available for plot_distribution")

  ratio_object <- simPower(
    n = 12, distribution = "lnorm", ctype = "ROM",
    mu_list = list(T = c(y1 = 1.05), R = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = .8)),
    list_uequi.tol = list(T_vs_R = c(y1 = 1.25)),
    nsim = 10, ncores = 1, seed = 2, keep_sim_data = TRUE
  )
  ratio_plot <- plot_distribution(ratio_object, estimand = "t_value",
                                  type = "histogram", n_trials = 5,
                                  seed = 2026)
  expect_s3_class(ratio_plot, "ggplot")
  expect_true(all(is.finite(ratio_plot$data$value)))
})

test_that("distribution plots use t-adjustment alpha", {
  skip_if_not_installed("ggplot2")
  object <- simPower(
    n = 12, distribution = "norm", ctype = "DOM",
    mu_list = list(T = c(y1 = 1.1, y2 = 1.1), R = c(y1 = 1, y2 = 1)),
    sigma_list = list(T = c(y1 = .2, y2 = .2), R = c(y1 = .2, y2 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.5, y2 = -.5)),
    list_uequi.tol = list(T_vs_R = c(y1 = .5, y2 = .5)),
    k = 1, adjust = "t", nsim = 20, ncores = 1, seed = 3,
    keep_sim_data = TRUE
  )
  t_plot <- plot_distribution(object, estimand = "t_value", n_trials = 10,
                              seed = 2026)
  expect_s3_class(t_plot, "ggplot")
  expect_gt(min(t_plot$data$critical[t_plot$data$boundary == "lower"]), 1.96)
})

test_that("overall diagnostic plots show the aggregate comparator result", {
  skip_if_not_installed("ggplot2")
  multi_object <- simPower(
    n = 12, distribution = "norm", ctype = "DOM",
    mu_list = list(T1 = c(y1 = 1.1), T2 = c(y1 = 1.05), R = c(y1 = 1)),
    sigma_list = list(T1 = c(y1 = .2), T2 = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T1_vs_R = c("T1", "R"),
                           T2_vs_R = c("T2", "R")),
    list_lequi.tol = list(T1_vs_R = c(y1 = -.5), T2_vs_R = c(y1 = -.5)),
    list_uequi.tol = list(T1_vs_R = c(y1 = .5), T2_vs_R = c(y1 = .5)),
    nsim = 10, ncores = 1, seed = 3
  )
  stability <- plot_stability(multi_object, overall = TRUE)
  mc_error <- plot_mc_error(multi_object, overall = TRUE)
  expect_equal(unique(as.character(stability$data$Comparator)), "All comparators")
  expect_equal(unique(as.character(mc_error$data$Comparator)), "All comparators")
  expect_s3_class(stability, "ggplot")
  expect_s3_class(mc_error, "ggplot")
})
