test_that("single-comparison count methods support Poisson and negative binomial models", {
  for (model in c("poisson", "negative-binomial")) {
    result <- SimTOST:::power_count(
      n_per_arm = 25, rate_test = 5.1, rate_reference = 5,
      exposure = 2, model = model, dispersion = 0.1,
      nsim = 30, seed = 11, ncores = 1
    )
    expect_s3_class(result, "countpower")
    expect_true(is.finite(result$power))
    expect_gte(result$power, 0)
    expect_lte(result$power, 1)
    expect_equal(result$n_total, 50)
    expect_named(confint(result), c("Lower", "Upper"))
    expect_true(is.data.frame(summary(result)))
  }
})

test_that("count sample-size methods return canonical count classes", {
  result <- SimTOST:::sampleSize_count(
    power = 0.8, rate_test = 5.1, rate_reference = 5,
    exposure = 2, model = "poisson", lower = 10, upper = 40,
    nsim = 20, seed = 12, ncores = 1
  )
  expect_s3_class(result, "countss")
  expect_true(result$n_per_arm >= 10)
  expect_true(result$n_per_arm <= 40)
  expect_true(is.data.frame(result$table.iter))
  expect_true(is.data.frame(result$table.test))
  expect_named(result$table.iter,
               c("n_iter", "n_total", "power", "power_LCI", "power_UCI"))
  expect_named(result$table.test,
               c("n_iter", "n_total", "totalyComp:comparison_1",
                 "endpoint_1Comp:comparison_1", "totaly"))
  expect_true(result$n_per_arm %in% result$table.iter$n_iter)
  expect_equal(nrow(result$table.test), nrow(result$table.iter) * 20)
  count_plot <- ggplot2::ggplot_build(plot(result))
  expect_s3_class(plot(result), "ggplot")
  expect_equal(nrow(count_plot$data[[2]]), nrow(result$table.iter) * 2)
  expect_true(all(count_plot$data[[3]]$linewidth >= 1))
  expect_true(is.data.frame(summary(result)))
  expect_named(confint(result), c("Lower", "Upper"))
})

test_that("count sample-size history handles a lower-bound solution and invalid history controls", {
  result <- SimTOST:::sampleSize_count(
    power = .1, rate_test = 5, rate_reference = 5,
    margin_lower = .1, margin_upper = 10, lower = 2, upper = 10,
    nsim = 5, seed = 121, ncores = 1
  )
  expect_equal(nrow(result$table.iter), 1L)
  expect_equal(unique(result$table.test$n_iter), 2)

  unified <- sampleSize(
    power = .1, distribution = "poisson",
    rate_list = list(TEST = c(y1 = 5), REF = c(y1 = 5)),
    list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
    list_lequi.tol = list(TEST_vs_REF = .1),
    list_uequi.tol = list(TEST_vs_REF = 10),
    lower = 2, upper = 2, nsim = 5, seed = 122, ncores = 1
  )
  expect_s3_class(unified, "simss")
  expect_s3_class(unified, "countss")
  expect_true(is.data.frame(unified$table.iter))
  expect_true(is.data.frame(unified$table.test))

  expect_error(
    SimTOST:::find_count_sample_size(
      power_fun = function(n) list(power = .5), target_power = .8,
      lower = 2, upper = 4, return_history = NA
    ),
    "return_history"
  )
})

test_that("unified count sample-size results match the continuous result contract", {
  result <- suppressWarnings(sampleSize(
    power = .1, distribution = "poisson",
    rate_list = list(TEST = c(y1 = 5), REF = c(y1 = 5)),
    list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
    list_y_comparator = list(TEST_vs_REF = "y1"),
    list_lequi.tol = list(TEST_vs_REF = .1),
    list_uequi.tol = list(TEST_vs_REF = 10),
    lower = 2, upper = 2, nsim = 5, seed = 123, ncores = 1
  ))

  expect_s3_class(result, "simss")
  expect_s3_class(result, "countss")
  expect_true(all(c("response", "table.iter", "table.test", "param.u",
                    "param", "param.d") %in% names(result)))
  expect_named(result$response,
               c("n_iter", "n_drop", "n_TEST", "n_REF", "n_total",
                 "power", "power_LCI", "power_UCI"))
  expect_equal(result$response$n_total, result$n_total)
  expect_equal(result$param.d$k, 1L)
  expect_identical(result$param$list_y_comparator[["TEST_vs_REF"]], "y1")
  expect_identical(result$parameters$.function, "sampleSize")

  printed <- capture.output(summary_result <- summary(result))
  expect_true(any(grepl("Sample Size Summary", printed, fixed = TRUE)))
  expect_true(any(grepl("Estimated Sample Size", printed, fixed = TRUE)))
  expect_true(any(grepl("event-rate ratio", printed, fixed = TRUE)))
  expect_named(summary_result, c("TEST", "REF", "Total"))
  expect_error(SimTOST:::summary.countss(list()), "class 'countss'")
})

test_that("joint count methods support multiple comparisons and endpoints", {
  rates <- list(TEST = c(y1 = 5, y2 = 6),
                REF1 = c(y1 = 5, y2 = 6),
                REF2 = c(y1 = 5, y2 = 6))
  comparisons <- list(TEST_vs_REF1 = c("TEST", "REF1"),
                      TEST_vs_REF2 = c("TEST", "REF2"))
  lower <- list(TEST_vs_REF1 = c(y1 = .8, y2 = .8),
                TEST_vs_REF2 = c(y1 = .8, y2 = .8))
  upper <- list(TEST_vs_REF1 = c(y1 = 1.25, y2 = 1.25),
                TEST_vs_REF2 = c(y1 = 1.25, y2 = 1.25))
  result <- SimTOST:::power_count_joint(
    n_per_arm = 20, rates = rates, comparisons = comparisons,
    list_margin_lower = lower, list_margin_upper = upper,
    endpoint_corr = diag(2), nsim = 20, seed = 13, ncores = 1
  )
  expect_s3_class(result, "countpower")
  expect_false(inherits(result, "countpower_joint"))
  expect_equal(result$n_comparisons, 2)
  expect_equal(result$n_endpoints, 2)
  expect_true(is.finite(result$power))
})

test_that("count methods reject invalid rates, margins, and controls", {
  expect_error(SimTOST:::power_count(
    20, rate_test = 0, rate_reference = 5, nsim = 10
  ), "Rates and exposure")
  expect_error(SimTOST:::power_count(
    20, rate_test = 5, rate_reference = 5,
    margin_lower = 1.25, margin_upper = .8, nsim = 10
  ), "Margins")
  expect_error(SimTOST:::power_count(
    20, rate_test = 5, rate_reference = 5, alpha = 1, nsim = 10
  ), "alpha")
  expect_error(SimTOST:::power_count(
    20, rate_test = 5, rate_reference = 5, nsim = 0
  ), "nsim")
})

test_that("count power uses selected endpoints and caps oversized k", {
  rates <- list(TEST = c(y1 = 5, y2 = 5, y3 = 5),
                REF = c(y1 = 5, y2 = 5, y3 = 5))
  comparison <- list(TEST_vs_REF = c("TEST", "REF"))
  lower <- list(TEST_vs_REF = c(y1 = .8, y2 = .8))
  upper <- list(TEST_vs_REF = c(y1 = 1.25, y2 = 1.25))

  fit <- NULL
  expect_warning(fit <- simPower(
    n = 20, distribution = "poisson", rate_list = rates,
    list_comparator = comparison,
    list_y_comparator = list(TEST_vs_REF = c("y1", "y2")),
    list_lequi.tol = lower, list_uequi.tol = upper,
    k = 3, adjust = "none", nsim = 20, seed = 14, ncores = 1
  ), "larger than the number")
  expect_equal(fit$n_endpoints, 2L)
  expect_equal(fit$k, 2)

  expect_error(simPower(
    n = 20, distribution = "poisson", rate_list = rates,
    list_comparator = comparison,
    list_y_comparator = list(TEST_vs_REF = "not_an_endpoint"),
    list_lequi.tol = list(TEST_vs_REF = .8),
    list_uequi.tol = list(TEST_vs_REF = 1.25),
    nsim = 5, ncores = 1
  ), "not available")

  inferred <- simPower(
    n = 20, distribution = "poisson", rate_list = rates,
    list_comparator = comparison,
    list_lequi.tol = list(TEST_vs_REF = c(y1 = .8, y2 = .8, y3 = .8)),
    list_uequi.tol = list(TEST_vs_REF = c(y1 = 1.25, y2 = 1.25, y3 = 1.25)),
    nsim = 5, ncores = 1, seed = 17
  )
  expect_equal(inferred$n_endpoints, 3L)
  expect_equal(inferred$k, 3)

  unnamed_type1 <- suppressWarnings(type1Error(
    null = "lower", joint = TRUE, n = 10, distribution = "poisson",
    rate_list = list(TEST = c(5, 5), REF = c(5, 5)),
    exposure = list(TEST = c(1, 1), REF = c(1, 1)),
    list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
    list_lequi.tol = list(TEST_vs_REF = c(.8, .8)),
    list_uequi.tol = list(TEST_vs_REF = c(1.25, 1.25)),
    nsim = 2, ncores = 1, seed = 18
  ))
  expect_s3_class(unnamed_type1, "type1error_joint")
  expect_setequal(unnamed_type1$table$Endpoint,
                  c("endpoint_1", "endpoint_2", "endpoint_1+endpoint_2"))
})

test_that("standalone count kernels cap oversized k and reject invalid k", {
  fit <- NULL
  expect_warning(fit <- SimTOST:::power_count(
    n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), k = 3,
    nsim = 10, seed = 15, ncores = 1
  ), "larger than the number")
  expect_equal(fit$k, 2)
  expect_error(SimTOST:::power_count(
    n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), k = c(1, 2), nsim = 5
  ), "positive integers")

  rates <- list(TEST = c(y1 = 5, y2 = 5), REF = c(y1 = 5, y2 = 5),
                ALT = c(y1 = 5, y2 = 5))
  comparisons <- list(TEST_vs_REF = c("TEST", "REF"),
                      TEST_vs_ALT = c("TEST", "ALT"))
  fit_joint <- NULL
  expect_warning(fit_joint <- SimTOST:::power_count_joint(
    n_per_arm = 20, rates = rates, comparisons = comparisons,
    list_margin_lower = list(TEST_vs_REF = c(y1 = .8, y2 = .8),
                             TEST_vs_ALT = c(y1 = .8, y2 = .8)),
    list_margin_upper = list(TEST_vs_REF = c(y1 = 1.25, y2 = 1.25),
                             TEST_vs_ALT = c(y1 = 1.25, y2 = 1.25)),
    k = 3, nsim = 10, seed = 16, ncores = 1
  ), "larger than the number")
  expect_equal(fit_joint$k, 2)

  expect_error(simPower(
    n = 20, distribution = "poisson",
    rate_list = list(TEST = c(y1 = 5, y2 = 5, y3 = 5),
                     REF = c(y1 = 5, y2 = 5, y3 = 5),
                     ALT = c(y1 = 5, y2 = 5, y3 = 5)),
    list_comparator = list(TEST_vs_REF = c("TEST", "REF"),
                           TEST_vs_ALT = c("TEST", "ALT")),
    list_y_comparator = list(TEST_vs_REF = c("y1", "y2"),
                             TEST_vs_ALT = c("y1", "y3")),
    list_lequi.tol = list(TEST_vs_REF = c(.8, .8),
                          TEST_vs_ALT = c(.8, .8)),
    list_uequi.tol = list(TEST_vs_REF = c(1.25, 1.25),
                          TEST_vs_ALT = c(1.25, 1.25)),
    nsim = 5, ncores = 1
  ), "same endpoints")
})

test_that("count Bonferroni warns but remains available when k equals m", {
  fit <- NULL
  expect_warning(fit <- SimTOST:::power_count(
    n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), margin_lower = .1,
    margin_upper = 10, k = 2,
    adjust = "bonferroni", nsim = 10, seed = 19, ncores = 1
  ), "not necessary")
  expect_equal(fit$k, 2)
  expect_identical(fit$adjust, "bonferroni")

  rates <- list(TEST = c(y1 = 5, y2 = 5), REF = c(y1 = 5, y2 = 5),
                ALT = c(y1 = 5, y2 = 5))
  comparisons <- list(TEST_vs_REF = c("TEST", "REF"),
                      TEST_vs_ALT = c("TEST", "ALT"))
  lower <- list(TEST_vs_REF = c(y1 = .1, y2 = .1),
                TEST_vs_ALT = c(y1 = .1, y2 = .1))
  upper <- list(TEST_vs_REF = c(y1 = 10, y2 = 10),
                TEST_vs_ALT = c(y1 = 10, y2 = 10))
  joint <- NULL
  expect_warning(joint <- SimTOST:::power_count_joint(
    n_per_arm = 20, rates = rates, comparisons = comparisons,
    list_margin_lower = lower, list_margin_upper = upper,
    k = 2, adjust = "bonferroni", nsim = 10, seed = 20, ncores = 1
  ), "not necessary")
  expect_equal(joint$k, 2)
})

test_that("count sample-size searches warn once for redundant Bonferroni", {
  fit <- NULL
  expect_warning(fit <- SimTOST:::sampleSize_count(
    power = .1, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), margin_lower = .1,
    margin_upper = 10, k = 2,
    adjust = "bonferroni", lower = 2, upper = 40,
    optimization_method = "fast", nsim = 5, seed = 21, ncores = 1
  ), "not necessary")
  expect_equal(fit$k, 2)

  rates <- list(TEST = c(y1 = 5, y2 = 5), REF = c(y1 = 5, y2 = 5),
                ALT = c(y1 = 5, y2 = 5))
  comparisons <- list(TEST_vs_REF = c("TEST", "REF"),
                      TEST_vs_ALT = c("TEST", "ALT"))
  lower <- list(TEST_vs_REF = c(y1 = .1, y2 = .1),
                TEST_vs_ALT = c(y1 = .1, y2 = .1))
  upper <- list(TEST_vs_REF = c(y1 = 10, y2 = 10),
                TEST_vs_ALT = c(y1 = 10, y2 = 10))
  joint <- NULL
  expect_warning(joint <- SimTOST:::sampleSize_count_joint(
    power = .1, rates = rates, comparisons = comparisons,
    list_margin_lower = lower, list_margin_upper = upper,
    k = 2, adjust = "bonferroni", lower = 2, upper = 40,
    optimization_method = "fast", nsim = 5, seed = 22, ncores = 1
  ), "not necessary")
  expect_equal(joint$k, 2)
  expect_true(is.data.frame(joint$table.iter))
  expect_true(is.data.frame(joint$table.test))
  expect_equal(nrow(joint$table.test), nrow(joint$table.iter) * 5)
  expect_true(all(c("totalyComp:TEST_vs_REF", "totalyComp:TEST_vs_ALT",
                    "totaly") %in% names(joint$table.test)))
  joint_plot <- ggplot2::ggplot_build(plot(joint, all = FALSE))
  expect_equal(nrow(joint_plot$layout$layout), 3L)
  joint_overall <- ggplot2::ggplot_build(plot(joint))
  expect_equal(nrow(joint_overall$layout$layout), 1L)
  expect_false("Comparator" %in% names(joint_overall$layout$layout))
  selected_plot <- ggplot2::ggplot_build(plot(joint, display = "TEST_vs_REF", all = FALSE))
  expect_equal(nrow(selected_plot$layout$layout), 1L)
  expect_error(plot(joint, display = "missing"), "Unknown comparator")
  expect_error(plot(joint, all = 1), "'all'")

  legacy <- joint
  legacy$table.test <- NULL
  expect_s3_class(plot(legacy), "ggplot")
})

test_that("count outcomes support sequential type_y for power and sample size", {
  fit <- SimTOST:::power_count(
    n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), margin_lower = .1,
    margin_upper = 10, k = 1, type_y = c(y1 = 1, y2 = 2),
    adjust = "seq", nsim = 10, seed = 23, ncores = 1
  )
  expect_s3_class(fit, "countpower")
  expect_identical(fit$adjust, "sequential")
  expect_identical(unname(fit$type_y), c(1L, 2L))
  expect_true(is.finite(fit$power))

  ss <- SimTOST:::sampleSize_count(
    power = .1, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), margin_lower = .1,
    margin_upper = 10, k = 1, type_y = c(y1 = 1, y2 = 2),
    adjust = "seq", lower = 2, upper = 20, nsim = 5, seed = 24,
    ncores = 1
  )
  expect_s3_class(ss, "countss")
  expect_equal(unname(ss$parameters$type_y), c(1L, 2L))
})

test_that("count outcomes support the t-adjustment", {
  expect_equal(SimTOST:::.normalize_count_adjustment("partial-conjunction"), "t")
  expect_equal(
    SimTOST:::.count_endpoint_alpha(.05, m = 2, k = 1, adjust = "t"),
    c(.025, .025)
  )
  fit <- SimTOST:::power_count(
    n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), margin_lower = .1,
    margin_upper = 10, k = 1, adjust = "t", nsim = 10, seed = 30,
    ncores = 1
  )
  expect_s3_class(fit, "countpower")
  expect_identical(fit$adjust, "t")
  expect_true(is.finite(fit$power))
  expect_warning(
    SimTOST:::power_count(
      n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
      rate_reference = c(y1 = 5, y2 = 5), k = 2, adjust = "t",
      nsim = 5, seed = 31, ncores = 1
    ),
    "has no effect"
  )
})

test_that("joint and unified count analyses support sequential type_y", {
  rates <- list(TEST = c(y1 = 5, y2 = 5),
                REF1 = c(y1 = 5, y2 = 5),
                REF2 = c(y1 = 5, y2 = 5))
  comparisons <- list(TEST_vs_REF1 = c("TEST", "REF1"),
                      TEST_vs_REF2 = c("TEST", "REF2"))
  lower <- list(TEST_vs_REF1 = c(y1 = .1, y2 = .1),
                TEST_vs_REF2 = c(y1 = .1, y2 = .1))
  upper <- list(TEST_vs_REF1 = c(y1 = 10, y2 = 10),
                TEST_vs_REF2 = c(y1 = 10, y2 = 10))
  fit <- SimTOST:::power_count_joint(
    n_per_arm = 20, rates = rates, comparisons = comparisons,
    list_margin_lower = lower, list_margin_upper = upper, k = 1,
    type_y = c(y1 = 1, y2 = 2), adjust = "seq", nsim = 10,
    seed = 25, ncores = 1
  )
  expect_identical(fit$adjust, "sequential")
  expect_identical(unname(fit$type_y), c(1L, 2L))

  unified <- simPower(
    n = 20, distribution = "poisson", rate_list = rates,
    list_comparator = comparisons,
    list_y_comparator = list(TEST_vs_REF1 = c("y1", "y2"),
                              TEST_vs_REF2 = c("y1", "y2")),
    list_lequi.tol = lower, list_uequi.tol = upper, k = 1,
    type_y = c(y1 = 1, y2 = 2), adjust = "seq", nsim = 10,
    seed = 26, ncores = 1
  )
  expect_s3_class(unified, "simpower")
  expect_identical(unified$parameters$adjust, "seq")
  expect_equal(unname(unified$parameters$type_y), c(1L, 2L))
})

test_that("count sequential hierarchy validates type_y and warns when absent", {
  expect_warning(
    SimTOST:::power_count(
      n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
      rate_reference = c(y1 = 5, y2 = 5), k = 1, adjust = "seq",
      nsim = 5, seed = 27, ncores = 1
    ),
    "requires 'type_y'"
  )
  expect_error(
    SimTOST:::power_count(
      n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
      rate_reference = c(y1 = 5, y2 = 5), k = 1,
      type_y = c(y1 = 1, y2 = 3), adjust = "seq", nsim = 5
    ),
    "only integer values 1"
  )
  expect_warning(
    SimTOST:::power_count(
      n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
      rate_reference = c(y1 = 5, y2 = 5), k = 1,
      type_y = c(y1 = 1, y2 = 2), adjust = "bonferroni", nsim = 5,
      seed = 28, ncores = 1
    ),
    "used only with adjust = 'seq'"
  )
})

test_that("type1Error can reuse a unified sequential count object", {
  fit <- simPower(
    n = 20, distribution = "poisson",
    rate_list = list(TEST = c(y1 = 5, y2 = 5), REF = c(y1 = 5, y2 = 5)),
    list_comparator = list(TEST_REF = c("TEST", "REF")),
    list_y_comparator = list(TEST_REF = c("y1", "y2")),
    list_lequi.tol = list(TEST_REF = c(y1 = .1, y2 = .1)),
    list_uequi.tol = list(TEST_REF = c(y1 = 10, y2 = 10)),
    k = 1, type_y = c(y1 = 1, y2 = 2), adjust = "seq",
    nsim = 5, seed = 29, ncores = 1
  )
  result <- type1Error(
    x = fit, null = "lower", joint = TRUE, nsim = 5, seed = 30,
    ncores = 1
  )
  expect_s3_class(result, "type1error_joint")
})
