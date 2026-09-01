# Re-run a SimTOST calculation while changing only explicitly supplied inputs.
.update_simtost_object <- function(object, fun, base_args, dots, evaluate) {
  if (!isTRUE(evaluate)) {
    return(as.call(c(list(as.name(fun)), utils::modifyList(base_args, dots))))
  }
  do.call(fun, utils::modifyList(base_args, dots))
}

.continuous_sample_args_from_object <- function(object) {
  result <- object
  if (inherits(object, "simpower")) result <- object$result
  if (is.null(result$param) || is.null(result$param.d))
    stop("The object does not retain the parameters needed for update().")
  p <- result$param
  d <- result$param.d
  mu_list <- lapply(seq_along(p$mu), function(i) {
    value <- as.numeric(p$mu[[i]])
    names(value) <- p$ynames_list[[i]]
    value
  })
  names(mu_list) <- names(p$mu)
  table_iter <- result$table.iter
  lower <- if (!is.null(d$lower)) d$lower else
    min(table_iter$n_iter, na.rm = TRUE)
  upper <- if (!is.null(d$upper)) d$upper else
    max(table_iter$n_iter, na.rm = TRUE)
  list(
    distribution = d$distribution,
    mu_list = mu_list,
    varcov_list = p$varcov,
    sigmaB = p$sigmaB,
    Eper = p$Eper,
    Eco = p$Eco,
    TAR = unlist(p$TAR_list),
    arm_names = p$arm_names,
    ynames_list = p$ynames_list,
    type_y = p$type_y,
    list_comparator = p$list_comparator,
    list_y_comparator = p$list_y_comparator,
    power = d$power,
    alpha = d$alpha,
    list_lequi.tol = d$list_lequi.tol,
    list_uequi.tol = d$list_uequi.tol,
    dtype = d$dtype,
    ctype = d$ctype,
    vareq = d$vareq,
    k = d$k,
    adjust = d$adjust,
    dropout = d$dropout,
    nsim = d$nsim,
    seed = if (!is.null(d$seed)) d$seed else 1234,
    ncores = if (!is.null(d$ncores)) d$ncores else 1,
    optimization_method = if (!is.null(d$optimization_method))
      d$optimization_method else "step-by-step",
    lower = lower,
    upper = upper,
    step.power = if (!is.null(d$step.power)) d$step.power else 6,
    step.up = if (!is.null(d$step.up)) d$step.up else TRUE,
    pos.side = if (!is.null(d$pos.side)) d$pos.side else FALSE,
    maxiter = if (!is.null(d$maxiter)) d$maxiter else 1000,
    verbose = if (!is.null(d$verbose)) d$verbose else FALSE
  )
}

.simpower_args_from_object <- function(object) {
  if (!is.null(object$parameters)) return(object$parameters)
  result <- object$result
  if (is.null(result))
    stop("The object does not retain the parameters needed for update().")
  args <- .continuous_sample_args_from_object(result)
  args$n <- object$n[[1L]]
  args$power <- NULL
  args$lower <- NULL
  args$upper <- NULL
  args$optimization_method <- NULL
  args$step.power <- NULL
  args$step.up <- NULL
  args$pos.side <- NULL
  args$maxiter <- NULL
  args$verbose <- NULL
  args$keep_sim_data <- !is.null(object$sim_data)
  args
}

.count_args_from_object <- function(object, fun) {
  if (!is.null(object$parameters)) {
    args <- object$parameters
    args$.function <- NULL
    return(args)
  }
  stop("The count object does not retain the original parameters needed for update().")
}

#' Update a SimTOST sample-size calculation
#'
#' Re-runs the calculation represented by `object`, replacing only arguments
#' supplied in `...`. This is useful for changing, for example, the target
#' power or number of simulations without repeating the complete original call.
#' @param object A `simss` object, including count-outcome sample-size objects.
#' @param ... Arguments to replace in the original calculation.
#' @param evaluate If `FALSE`, return the updated call instead of evaluating it.
#' @return A newly calculated object of the same result family as `object`, or
#'   an unevaluated call when `evaluate = FALSE`.
#' @export
#' @method update simss
update.simss <- function(object, ..., evaluate = TRUE) {
  if (!inherits(object, "simss")) stop("'object' must be a simss object.")
  if (inherits(object, "countss")) {
    updated <- update.countss(object, ..., evaluate = evaluate)
    if (!isTRUE(evaluate)) return(updated)
    class(updated) <- c("simss", class(updated))
    return(updated)
  }
  .update_simtost_object(
    object, "sampleSize", .continuous_sample_args_from_object(object),
    list(...), evaluate
  )
}

#' Update a SimTOST fixed-sample-size power calculation
#'
#' @param object A `simpower` object, including unified count-outcome results.
#' @param ... Arguments to replace in the original calculation.
#' @param evaluate If `FALSE`, return the updated call instead of evaluating it.
#' @return A newly calculated `simpower` object, or an unevaluated call.
#' @export
#' @method update simpower
update.simpower <- function(object, ..., evaluate = TRUE) {
  if (!inherits(object, "simpower")) stop("'object' must be a simpower object.")
  .update_simtost_object(
    object, "simPower", .simpower_args_from_object(object),
    list(...), evaluate
  )
}

#' Update a standalone count power calculation
#'
#' @param object A standalone `countpower` object returned by a count power
#'   function.
#' @param ... Arguments to replace in the original calculation.
#' @param evaluate If `FALSE`, return the updated call instead of evaluating it.
#' @return A newly calculated count power object, or an unevaluated call.
#' @export
#' @method update countpower
update.countpower <- function(object, ..., evaluate = TRUE) {
  if (!inherits(object, "countpower"))
    stop("'object' must be a countpower object.")
  args <- .count_args_from_object(object, object$parameters$.function)
  fun <- object$parameters$.function
  if (is.null(fun)) fun <- if (!is.null(object$n_comparisons))
    "power_count_joint" else "power_count"
  .update_simtost_object(object, fun, args, list(...), evaluate)
}

#' Update a standalone count sample-size calculation
#'
#' @param object A standalone `countss` object returned by a count sample-size
#'   function.
#' @param ... Arguments to replace in the original calculation.
#' @param evaluate If `FALSE`, return the updated call instead of evaluating it.
#' @return A newly calculated count sample-size object, or an unevaluated call.
#' @export
#' @method update countss
update.countss <- function(object, ..., evaluate = TRUE) {
  if (!inherits(object, "countss"))
    stop("'object' must be a countss object.")
  fun <- object$parameters$.function
  if (is.null(fun)) fun <- if (!is.null(object$parameters$comparisons))
    "sampleSize_count_joint" else "sampleSize_count"
  args <- .count_args_from_object(object, fun)
  .update_simtost_object(object, fun, args, list(...), evaluate)
}
