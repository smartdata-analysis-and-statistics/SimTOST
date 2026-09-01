#' Run a simulation for a supported design and test
#'
#' @description A single public dispatcher for the design/test combinations
#' implemented by SimTOST. The compiled design-specific routines remain
#' internal implementation details.
#' @param design Trial design: `"parallel"` or `"2x2"`.
#' @param ctype Test type: `"DOM"` or `"ROM"`.
#' @param ... Arguments passed to the corresponding internal simulation routine.
#' @return The result returned by the selected simulation routine.
#' @export
#' @examples
#' # The dispatcher is primarily used by the package internals. For example:
#' # run_simulations(design = "parallel", ctype = "DOM", ...)
run_simulations <- function(design = c("parallel", "2x2"),
                            ctype = c("DOM", "ROM"), ...) {
  design <- match.arg(design)
  ctype <- match.arg(ctype)
  fun <- switch(paste(design, ctype, sep = "_"),
    parallel_DOM = run_simulations_par_dom,
    parallel_ROM = run_simulations_par_rom,
    `2x2_DOM` = run_simulations_2x2_dom,
    `2x2_ROM` = run_simulations_2x2_rom
  )
  fun(...)
}
