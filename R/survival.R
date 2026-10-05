# Survival analysis helpers ----------------------------------------------

#' Fit Kaplan-Meier survival curves
#'
#' Fits `survfit(Surv(time, event) ~ group, data = ...)` and returns the
#' fitted object for use with [plot_km_curves()] or base plot methods.
#'
#' @importFrom survival Surv survfit coxph
#'
#' @param survival_data Data frame with time, event and grouping columns.
#' @param time Character. Column name of the follow-up time.
#' @param event Character. Column name of the event indicator
#'   (0/1 or FALSE/TRUE).
#' @param group Character. Column name of the grouping variable.
#'
#' @return A `survfit` object.
#'
#' @export
#' @examples
#' fit <- fit_km(example_survival, time = "OS", event = "status",
#'               group = "cluster")
#' fit
fit_km <- function(survival_data, time, event, group) {
  .check_survival_columns(survival_data, time, event, group)
  .fit_survival("survfit", survival_data, time, event, group)
}

#' Fit a Cox proportional hazards model on a grouping variable
#'
#' Fits `coxph(Surv(time, event) ~ group, data = ...)` and returns the model
#' object; call [summary()] on the result to obtain hazard ratios and
#' p-values.
#'
#' @param survival_data Data frame with time, event and grouping columns.
#' @param time Character. Column name of the follow-up time.
#' @param event Character. Column name of the event indicator.
#' @param group Character. Column name of the grouping variable.
#'
#' @return A `coxph` object.
#'
#' @export
#' @examples
#' fit <- fit_cox(example_survival, time = "OS", event = "status",
#'                group = "cluster")
#' summary(fit)
fit_cox <- function(survival_data, time, event, group) {
  .check_survival_columns(survival_data, time, event, group)
  .fit_survival("coxph", survival_data, time, event, group)
}

#' Plot Kaplan-Meier curves with risk table
#'
#' Wraps [survminer::ggsurvplot()] with sensible defaults (log-rank p-value,
#' risk table, NPG palette) as used in the original pipeline.
#'
#' @param fit A `survfit` object from [fit_km()].
#' @param survival_data Data frame used to fit the model.
#' @param ... Further arguments passed to [survminer::ggsurvplot()].
#'
#' @return A `ggsurvplot` object (print it to draw).
#'
#' @export
#' @examples
#' fit <- fit_km(example_survival, time = "OS", event = "status",
#'               group = "cluster")
#' plot_km_curves(fit, example_survival)
plot_km_curves <- function(fit, survival_data, ...) {
  survminer::ggsurvplot(
    fit,
    data = survival_data,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    palette = ggsci::pal_npg()(10),
    ggtheme = ggplot2::theme_classic(),
    legend = "right",
    xlab = "Time (Months)",
    ylab = "Survival Probability",
    ...
  )
}

.check_survival_columns <- function(survival_data, time, event, group) {
  needed <- c(time, event, group)
  missing <- setdiff(needed, colnames(survival_data))
  if (length(missing) > 0) {
    stop("`survival_data` is missing column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Fits `survival::<fun>(Surv(time, event) ~ group, data = ...)` keeping the
# formula inline in the call (required by survminer::ggsurvplot()) while
# making `Surv` resolvable from any calling environment.
.fit_survival <- function(fun, survival_data, time, event, group) {
  env <- new.env(parent = asNamespace("survival"))
  env$survival_data <- survival_data

  expr <- parse(text = sprintf(
    "survival::%s(Surv(%s, %s) ~ %s, data = survival_data)",
    fun, time, event, group
  ))

  eval(expr, env)
}
