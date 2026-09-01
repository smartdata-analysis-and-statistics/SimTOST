test_that("type1Error places one endpoint on each least-favorable boundary", {
  result <- simPower(
    n = 30,
    distribution = "normal",
    mu_list = list(R = c(y1 = 1, y2 = 1), T = c(y1 = 1, y2 = 1)),
    sigma_list = list(R = c(y1 = 0.2, y2 = 0.2),
                      T = c(y1 = 0.2, y2 = 0.2)),
    list_comparator = list(R_vs_T = c("R", "T")),
    list_lequi.tol = list(R_vs_T = c(y1 = -0.2, y2 = -0.2)),
    list_uequi.tol = list(R_vs_T = c(y1 = 0.2, y2 = 0.2)),
    ctype = "DOM", nsim = 30, ncores = 1, seed = 123
  )

  boundary <- type1Error(x = result, null = "both", nsim = 30, seed = 123)

  expect_s3_class(boundary, "type1error_set")
  expect_identical(boundary$lower$boundary_endpoint, "y1")
  expect_identical(boundary$upper$boundary_endpoint, "y1")
  expect_equal(boundary$lower$boundary_estimands$Estimand[1], -0.2)
  expect_equal(boundary$upper$boundary_estimands$Estimand[1], 0.2)
  expect_true(is.finite(boundary$lower$type1_error))
  expect_true(is.finite(boundary$upper$type1_error))
  expect_true(is.finite(attr(boundary, "global_upper")))
})

test_that("type1Error accepts an explicit endpoint and rejects unknown endpoints", {
  result <- simPower(
    n = 20,
    distribution = "normal",
    mu_list = list(R = c(y1 = 1, y2 = 1), T = c(y1 = 1, y2 = 1)),
    sigma_list = list(R = c(y1 = 0.2, y2 = 0.2),
                      T = c(y1 = 0.2, y2 = 0.2)),
    list_comparator = list(R_vs_T = c("R", "T")),
    list_lequi.tol = list(R_vs_T = c(y1 = -0.2, y2 = -0.2)),
    list_uequi.tol = list(R_vs_T = c(y1 = 0.2, y2 = 0.2)),
    ctype = "DOM", nsim = 20, ncores = 1, seed = 456
  )

  selected <- type1Error(x = result, null = "lower", comparator = "R_vs_T",
                         endpoint = "y2", nsim = 20, seed = 456)
  expect_identical(selected$boundary_comparator, "R_vs_T")
  expect_identical(selected$boundary_endpoint, "y2")
  expect_equal(selected$boundary_estimands$Estimand[
    selected$boundary_estimands$Endpoint == "y2"], -0.2)

  expect_error(
    type1Error(x = result, null = "lower", endpoint = "unknown",
               nsim = 10, seed = 456),
    "selected comparator"
  )
  expect_error(
    type1Error(x = result, null = "lower", comparator = "unknown",
               nsim = 10, seed = 456),
    "comparator"
  )
})

test_that("type1Error joint mode evaluates scenarios and plots results", {
  result <- simPower(
    n = 20,
    distribution = "normal",
    mu_list = list(R = c(y1 = 1, y2 = 1), T = c(y1 = 1, y2 = 1)),
    sigma_list = list(R = c(y1 = 0.2, y2 = 0.2),
                      T = c(y1 = 0.2, y2 = 0.2)),
    list_comparator = list(R_vs_T = c("R", "T")),
    list_lequi.tol = list(R_vs_T = c(y1 = -0.2, y2 = -0.2)),
    list_uequi.tol = list(R_vs_T = c(y1 = 0.2, y2 = 0.2)),
    ctype = "DOM", nsim = 20, ncores = 1, seed = 789
  )

  joint <- type1Error(x = result, null = "both", joint = TRUE,
                      nsim = 20, seed = 789)

  expect_s3_class(joint, "type1error_joint")
  expect_equal(nrow(joint$table), 8L)
  expect_setequal(unique(joint$table$Boundary),
                  c("lower", "upper", "lower/upper", "upper/lower"))
  expect_setequal(unique(joint$table$Endpoint),
                  c("y1", "y2", "y1+y2"))
  expect_false(any(c("Component_Type1_Error", "Component_Lower", "Component_Upper") %in% names(joint$table)))
  joint_plot <- plot(joint)
  expect_s3_class(joint_plot, "ggplot")
  built <- ggplot2::ggplot_build(joint_plot)
  expect_false(any(vapply(built$data, function(layer)
    "linetype" %in% names(layer) && any(layer$linetype == "dotted"),
    logical(1))))
  expect_s3_class(joint_plot$theme$strip.text, "element_blank")
  line_idx <- which(vapply(built$data, function(layer)
    all(c("ymin", "ymax", "linewidth") %in% names(layer)), logical(1)))[[1L]]
  point_idx <- which(vapply(built$data, function(layer)
    all(c("x", "y", "size") %in% names(layer)) &&
      !"ymin" %in% names(layer), logical(1)))[[1L]]
  expect_true(max(built$data[[line_idx]]$linewidth) >
                min(built$data[[line_idx]]$linewidth))
  expect_true(max(built$data[[point_idx]]$size) >
                min(built$data[[point_idx]]$size))
  label_layers <- vapply(built$data,
                         function(layer) "label" %in% names(layer),
                         logical(1))
  expect_true(any(label_layers))
  label_layer <- built$data[[which(label_layers)[[1L]]]]
  label_row <- which(nzchar(label_layer$label))[[1L]]
  point_row <- which.max(built$data[[point_idx]]$y)
  expect_identical(label_layer$label[[label_row]],
                   sprintf("%.1f%%", 100 * max(joint$table$Type1_Error)))
  expect_equal(label_layer$x[[label_row]], built$data[[point_idx]]$x[[point_row]])
  expect_gt(label_layer$y[[label_row]], built$data[[point_idx]]$y[[point_row]])
  expect_false(grepl("Component_", paste(capture.output(print(joint)), collapse = "\n")))
  expect_equal(nrow(joint$global), 1L)
  expect_true(is.finite(joint$global_upper))
  expect_true(all(joint$table$Simultaneous_Upper >= joint$table$Type1_Error))
  expect_error(type1Error(x = result, joint = "yes"), "joint")
  expect_error(type1Error(x = result, joint = TRUE, conf.level = 1),
               "conf.level")
  expect_error(plot.type1error_joint(list()), "type1error_joint")
})

test_that("single-boundary plots show joint and global errors", {
  result <- simPower(
    n = 10, distribution = "norm", ctype = "DOM",
    mu_list = list(R = c(y1 = 1), T = c(y1 = 1)),
    sigma_list = list(R = c(y1 = .2), T = c(y1 = .2)),
    list_comparator = list(R_vs_T = c("R", "T")),
    list_lequi.tol = list(R_vs_T = c(y1 = -.2)),
    list_uequi.tol = list(R_vs_T = c(y1 = .2)),
    nsim = 10, ncores = 1, seed = 321
  )
  one <- type1Error(x = result, null = "lower", nsim = 10, seed = 321)
  expect_false(grepl("boundary-component", paste(capture.output(print(one)), collapse = "\n")))
  one_plot <- plot(one)
  expect_s3_class(one_plot, "ggplot")
  one_built <- ggplot2::ggplot_build(one_plot)
  point_idx <- which(vapply(one_built$data, function(layer)
    all(c("x", "y", "size", "colour") %in% names(layer)) &&
      !"ymin" %in% names(layer), logical(1)))[[1L]]
  expect_identical(one_built$data[[point_idx]]$colour[[1L]],
                   .type1_boundary_colors[["lower"]])
  label_layers <- vapply(one_built$data,
                         function(layer) "label" %in% names(layer),
                         logical(1))
  expect_true(any(label_layers))
  label_layer <- one_built$data[[which(label_layers)[[1L]]]]
  expect_identical(label_layer$label[[1L]],
                   sprintf("%.1f%%", 100 * one$type1_error))
})

test_that("joint Type I error produces a global result across comparisons", {
  result <- simPower(
    n = 15,
    distribution = "norm", ctype = "DOM",
    mu_list = list(T1 = c(y1 = 1), T2 = c(y1 = 1), R = c(y1 = 1)),
    sigma_list = list(T1 = c(y1 = .2), T2 = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T1_vs_R = c("T1", "R"),
                           T2_vs_R = c("T2", "R")),
    list_lequi.tol = list(T1_vs_R = c(y1 = -.2), T2_vs_R = c(y1 = -.2)),
    list_uequi.tol = list(T1_vs_R = c(y1 = .2), T2_vs_R = c(y1 = .2)),
    nsim = 10, ncores = 1, seed = 123
  )
  joint <- type1Error(x = result, null = "both", joint = TRUE,
                      nsim = 10, seed = 123)

  expect_equal(nrow(joint$table), 4L)
  expect_equal(nrow(joint$global), 1L)
  expect_true(joint$global$Type1_Error[[1L]] >=
                max(joint$table$Type1_Error) - .Machine$double.eps)
  expect_s3_class(plot(joint), "ggplot")
})

test_that("log-normal ROM boundaries are placed on the analysis scale", {
  boundary <- type1Error(
    null = "lower", n = 20, distribution = "lnorm", ctype = "ROM",
    mu_list = list(T = c(y1 = 0.84), R = c(y1 = 1.05)),
    sigma_list = list(T = c(y1 = 0.20), R = c(y1 = 0.20)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = "y1"),
    list_lequi.tol = list(T_vs_R = c(y1 = 0.80)),
    list_uequi.tol = list(T_vs_R = c(y1 = 1.25)),
    k = 1, nsim = 20, seed = 123, ncores = 1
  )

  # The arithmetic mean required to put the log-scale ROM exactly at 0.80
  # is not 0.84 when the arm means differ but the SDs are equal.
  expect_equal(boundary$result$param$mu$T[[1L]], 0.8478135, tolerance = 1e-6)
  expect_equal(boundary$boundary_estimands$Estimand[[1L]], 0.80,
               tolerance = 1e-12)

  upper <- type1Error(
    null = "upper", n = 20, distribution = "lnorm", ctype = "ROM",
    mu_list = list(T = c(y1 = 1), R = c(y1 = 1.05)),
    sigma_list = list(T = c(y1 = 0.20), R = c(y1 = 0.20)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = "y1"),
    list_lequi.tol = list(T_vs_R = c(y1 = 0.80)),
    list_uequi.tol = list(T_vs_R = c(y1 = 1.25)),
    k = 1, nsim = 20, seed = 123, ncores = 1
  )
  expect_equal(upper$result$param$mu$T[[1L]], 1.304387, tolerance = 1e-6)
  expect_equal(upper$boundary_estimands$Estimand[[1L]], 1.25,
               tolerance = 1e-12)
})

test_that("joint Type1_Error uses k from the source simulation object", {
  result <- simPower(
    n = 20, distribution = "normal", ctype = "DOM",
    mu_list = list(T1 = c(y1 = 1, y2 = 1),
                   T2 = c(y1 = 1, y2 = 1),
                   R = c(y1 = 1, y2 = 1)),
    sigma_list = list(T1 = c(y1 = .2, y2 = .2),
                      T2 = c(y1 = .2, y2 = .2),
                      R = c(y1 = .2, y2 = .2)),
    list_comparator = list(T1_vs_R = c("T1", "R"),
                           T2_vs_R = c("T2", "R")),
    list_y_comparator = list(T1_vs_R = c("y1", "y2"),
                             T2_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T1_vs_R = c(y1 = -.2, y2 = -.2),
                         T2_vs_R = c(y1 = -.2, y2 = -.2)),
    list_uequi.tol = list(T1_vs_R = c(y1 = .2, y2 = .2),
                         T2_vs_R = c(y1 = .2, y2 = .2)),
    k = c(2, 2), nsim = 20, ncores = 1, seed = 456
  )
  joint <- type1Error(x = result, null = "lower", joint = TRUE,
                      nsim = 20, seed = 456)

  for (i in seq_len(nrow(joint$table))) {
    key <- with(joint$table[i, ], paste(Boundary, Comparator, Endpoint,
                                        sep = "__"))
    scenario <- joint$results[[key]]
    expect_equal(scenario$result$param.d$k, c(2, 2))
    expect_equal(joint$table$Type1_Error[[i]], scenario$type1_error)
    expect_equal(joint$table$Lower[[i]], scenario$power_LCI)
    expect_equal(joint$table$Upper[[i]], scenario$power_UCI)
    expect_equal(scenario$mc_successes,
                 sum(scenario$result[["table.test"]][["totaly"]]))
  }
})

test_that("log-normal ROM boundary construction rejects missing variances", {
  expect_error(
    type1Error(
      null = "lower", n = 20, distribution = "lnorm", ctype = "ROM",
      mu_list = list(T = c(y1 = 0.84), R = c(y1 = 1.05)),
      list_comparator = list(T_vs_R = c("T", "R")),
      list_y_comparator = list(T_vs_R = "y1"),
      list_lequi.tol = list(T_vs_R = c(y1 = 0.80)),
      list_uequi.tol = list(T_vs_R = c(y1 = 1.25)),
      k = 1, nsim = 10, ncores = 1
    ),
    "variances"
  )
})

test_that("log-normal nonselected components use the calibrated midpoint", {
  joint <- type1Error(
    null = "lower", joint = TRUE, n = 10, distribution = "lnorm",
    ctype = "ROM",
    mu_list = list(T = c(y1 = 1, y2 = 1),
                   R = c(y1 = 1.05, y2 = 1.05)),
    sigma_list = list(T = c(y1 = .2, y2 = .3),
                      R = c(y1 = .2, y2 = .1)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = .8, y2 = .8)),
    list_uequi.tol = list(T_vs_R = c(y1 = 1.25, y2 = 1.25)),
    k = 2, nsim = 10, ncores = 1, seed = 321
  )
  scenario <- joint$results[["lower__T_vs_R__y1"]]
  target <- scenario$result$param$mu$T[1, "y2"]
  reference <- scenario$result$param$mu$R[1, "y2"]
  analysis_ratio <- exp(
    log(target) - .5 * log1p(.3^2 / target^2) -
      log(reference) + .5 * log1p(.1^2 / reference^2)
  )

  expect_equal(unname(analysis_ratio), sqrt(.8 * 1.25), tolerance = 1e-8)
  expect_gt(abs(target / reference - 1), 1e-3)
})

test_that("type1Error uses user-supplied lower and upper margins", {
  dom_lower <- type1Error(
    null = "lower", n = 10, distribution = "norm", ctype = "DOM",
    mu_list = list(T = c(y1 = 1), R = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = "y1"),
    list_lequi.tol = list(T_vs_R = c(y1 = -.4)),
    list_uequi.tol = list(T_vs_R = c(y1 = .6)),
    nsim = 10, ncores = 1, seed = 11
  )
  dom_upper <- type1Error(
    null = "upper", n = 10, distribution = "norm", ctype = "DOM",
    mu_list = list(T = c(y1 = 1), R = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = "y1"),
    list_lequi.tol = list(T_vs_R = c(y1 = -.4)),
    list_uequi.tol = list(T_vs_R = c(y1 = .6)),
    nsim = 10, ncores = 1, seed = 11
  )
  expect_equal(dom_lower$boundary_estimands$Estimand[[1L]], -.4)
  expect_equal(dom_upper$boundary_estimands$Estimand[[1L]], .6)

  count_lower <- type1Error(
    null = "lower", n = 10, distribution = "pois", ctype = "RR",
    rate_list = list(T = c(y1 = 2), R = c(y1 = 2)),
    exposure = list(T = c(y1 = 1), R = c(y1 = 1)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = "y1"),
    list_lequi.tol = list(T_vs_R = c(y1 = .6)),
    list_uequi.tol = list(T_vs_R = c(y1 = 1.4)),
    nsim = 10, ncores = 1, seed = 12
  )
  expect_equal(count_lower$boundary_estimands$Estimand[[1L]], .6)
})

test_that("joint type1Error enumerates k-aware composite null configurations", {
  make_source <- function(k) {
    suppressWarnings(simPower(
      n = 8, distribution = "norm", ctype = "DOM",
      mu_list = list(
        T = c(y1 = 1, y2 = 1, y3 = 1),
        R1 = c(y1 = 1, y2 = 1, y3 = 1),
        R2 = c(y1 = 1, y2 = 1, y3 = 1)
      ),
      sigma_list = list(
        T = c(y1 = .2, y2 = .2, y3 = .2),
        R1 = c(y1 = .2, y2 = .2, y3 = .2),
        R2 = c(y1 = .2, y2 = .2, y3 = .2)
      ),
      list_comparator = list(
        T_vs_R1 = c("T", "R1"),
        T_vs_R2 = c("T", "R2")
      ),
      list_y_comparator = list(
        T_vs_R1 = c("y1", "y2", "y3"),
        T_vs_R2 = c("y1", "y2", "y3")
      ),
      list_lequi.tol = list(
        T_vs_R1 = c(y1 = -.2, y2 = -.2, y3 = -.2),
        T_vs_R2 = c(y1 = -.2, y2 = -.2, y3 = -.2)
      ),
      list_uequi.tol = list(
        T_vs_R1 = c(y1 = .2, y2 = .2, y3 = .2),
        T_vs_R2 = c(y1 = .2, y2 = .2, y3 = .2)
      ),
      k = c(k, k), nsim = 2, ncores = 1, seed = 901
    ))
  }

  k3 <- suppressWarnings(type1Error(x = make_source(3), null = "both", joint = TRUE,
                                    nsim = 2, seed = 902))
  expect_equal(nrow(k3$table), 52L)
  expect_true(all(k3$table$NullCount %in% 1:3))
  expect_true(all(vapply(k3$results, function(z) {
    nrow(z$boundary_estimands) == z$boundary_null_count
  }, logical(1))))

  k2 <- suppressWarnings(type1Error(x = make_source(2), null = "both", joint = TRUE,
                                    nsim = 2, seed = 903))
  expect_equal(nrow(k2$table), 40L)
  expect_setequal(unique(k2$table$Endpoint),
                  c("y1+y2", "y1+y3", "y2+y3", "y1+y2+y3"))
  expect_true(any(grepl("/", k2$table$Boundary, fixed = TRUE)))
  expect_s3_class(plot(k2), "ggplot")
  minimal_plot <- ggplot2::ggplot_build(plot(k2))
  minimal_points <- which(vapply(
    minimal_plot$data,
    function(layer) all(c("x", "y", "size") %in% names(layer)) &&
      !"ymin" %in% names(layer),
    logical(1)
  ))[[1L]]
  expect_equal(
    nrow(minimal_plot$data[[minimal_points]]),
    sum(k2$table$NullCount == min(k2$table$NullCount))
  )
  all_plot <- ggplot2::ggplot_build(plot(k2, null_count = "all"))
  all_points <- which(vapply(
    all_plot$data,
    function(layer) all(c("x", "y", "size") %in% names(layer)) &&
      !"ymin" %in% names(layer),
    logical(1)
  ))[[1L]]
  expect_equal(nrow(all_plot$data[[all_points]]), nrow(k2$table))
  expect_error(plot(k2, null_count = "invalid"), "arg")
  expect_true(all(vapply(k2$results, function(z) {
    nrow(z$boundary_estimands) == z$boundary_null_count
  }, logical(1))))

  k1 <- suppressWarnings(type1Error(x = make_source(1), null = "both", joint = TRUE,
                                    nsim = 2, seed = 904))
  expect_equal(nrow(k1$table), 16L)
  expect_true(all(k1$table$Endpoint == "y1+y2+y3"))
  expect_true(all(vapply(k1$results, function(z) {
    nrow(z$boundary_estimands) == 3L
  }, logical(1))))

  four_source <- suppressWarnings(simPower(
    n = 8, distribution = "norm", ctype = "DOM",
    mu_list = list(T = setNames(rep(1, 4), paste0("y", 1:4)),
                   R = setNames(rep(1, 4), paste0("y", 1:4))),
    sigma_list = list(T = setNames(rep(.2, 4), paste0("y", 1:4)),
                      R = setNames(rep(.2, 4), paste0("y", 1:4))),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = paste0("y", 1:4)),
    list_lequi.tol = list(T_vs_R = setNames(rep(-.2, 4), paste0("y", 1:4))),
    list_uequi.tol = list(T_vs_R = setNames(rep(.2, 4), paste0("y", 1:4))),
    k = 4, nsim = 2, ncores = 1, seed = 905
  ))
  four <- suppressWarnings(type1Error(x = four_source, null = "both",
                                      joint = TRUE, nsim = 2, seed = 906))
  expect_equal(nrow(four$table), 80L)
  expect_true(all(vapply(four$results, function(z) {
    nrow(z$boundary_estimands) == z$boundary_null_count
  }, logical(1))))
  capped <- NULL
  expect_warning(
    capped <- type1Error(x = four_source, null = "both", joint = TRUE, k = 5,
                         nsim = 2, seed = 906),
    "larger than the number"
  )
  expect_equal(unname(capped$results[[1]]$result$param.d$k), 4L)
  expect_true(is.finite(capped$global_upper))
})

test_that("composite type1Error scenarios require joint mode", {
  expect_error(
    type1Error(
      null = "lower", joint = FALSE, n = 8, distribution = "norm",
      ctype = "DOM",
      mu_list = list(T = c(y1 = 1, y2 = 1, y3 = 1),
                     R = c(y1 = 1, y2 = 1, y3 = 1)),
      sigma_list = list(T = c(y1 = .2, y2 = .2, y3 = .2),
                        R = c(y1 = .2, y2 = .2, y3 = .2)),
      list_comparator = list(T_vs_R = c("T", "R")),
      list_y_comparator = list(T_vs_R = c("y1", "y2", "y3")),
      list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2, y3 = -.2)),
      list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2, y3 = .2)),
      k = 2, nsim = 2, ncores = 1
    ),
    "joint = TRUE"
  )
})

test_that("type1Error evaluates the t-adjustment", {
  source <- suppressWarnings(simPower(
    n = 8, distribution = "norm", ctype = "DOM",
    mu_list = list(T = c(y1 = 1, y2 = 1, y3 = 1),
                   R = c(y1 = 1, y2 = 1, y3 = 1)),
    sigma_list = list(T = c(y1 = .2, y2 = .2, y3 = .2),
                      R = c(y1 = .2, y2 = .2, y3 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2", "y3")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2, y3 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2, y3 = .2)),
    k = 2, adjust = "t", nsim = 2, ncores = 1, seed = 907
  ))
  result <- suppressWarnings(type1Error(
    x = source, null = "both", joint = TRUE, nsim = 2, seed = 908
  ))
  expect_equal(nrow(result$table), 20L)
  expect_true(all(vapply(result$results, function(z) {
    identical(z$result$param.d$adjust, "t")
  }, logical(1))))
  expect_true(is.finite(result$global_upper))
  expect_true(all(is.finite(result$table$Simultaneous_Upper)))
})
