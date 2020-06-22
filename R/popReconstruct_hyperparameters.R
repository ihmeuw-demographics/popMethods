#' @title Calculate Hyperparameters for the popReconstruct Model
#'
#' @description
#' Use elicited expert opinion about the observable marginal quantities in the
#' popReconstruct model to quantify uncertainty in the input parameters.
#'
#' Elicitation statement is "there is a \eqn{p_{\nu}}% probability that the true
#' values are within plus-or-minus \eqn{\eta_{\nu}} x 100 percent of the initial
#' point estimates. See references for more details.
#'
#' @param abs_deviations \[`list()`\]\cr
#'   Named list of absolute deviations \eqn{\eta_{\nu}} used in elicitation
#'   statement for each ccmpp/model component.
#' @param p \[`list()`\]\cr
#'   Named list of central \eqn{p_{\nu}}% probability interval used in
#'   elicitation statement for each ccmpp/model component.
#' @param alpha \[`list()`\]\cr
#'   Named list for alpha (shape) parameter of inverse gamma distribution for
#'   each ccmpp/model component.
#' @param survival \[`numeric()`\]\cr
#'   Initial age-sex specific survival proportion estimates. Must be non-NULL
#'   if 'survival' parameter is included in the other input arguments.
#'
#' @return Nested list of alpha (shape) and beta (scale) parameters for the
#'   variance of each model component.
#'
#' @inherit popReconstruct_fit references
#'
#' @family popReconstruct
#'
#' @examples
#' abs_deviations <- list(
#'   srb = 0.1,
#'   asfr = 0.1,
#'   population = 0.1,
#'   survival = 0.1,
#'   net_migration = 0.2
#' )
#' p <- list(
#'   srb = 0.9,
#'   asfr = 0.9,
#'   population = 0.9,
#'   survival = 0.9,
#'   net_migration = 0.9
#' )
#' alpha <- list(
#'   srb = 0.5,
#'   asfr = 0.5,
#'   population = 0.5,
#'   survival = 0.5,
#'   net_migration = 0.5
#' )
#' hyperparameters <- popReconstruct_hyperparameters(
#'   abs_deviations,
#'   p,
#'   alpha,
#'   survival = demCore::thailand_initial_estimates$survival$value
#' )
#'
#' @export
popReconstruct_hyperparameters <- function(abs_deviations,
                                           p,
                                           alpha,
                                           survival = NULL) {

  # Validate arguments ------------------------------------------------------

  possible_components <- c("srb", "asfr", "baseline", "population", "survival",
                           "net_migration", "immigration", "emigration")
  log_transform_components <- setdiff(possible_components,
                                      c("survival", "net_migration"))

  assertthat::assert_that(
    assertive::is_list(abs_deviations),
    assertive::is_list(p),
    assertive::is_list(alpha),
    all(sapply(abs_deviations, assertive::is_numeric)),
    all(sapply(p, assertive::is_numeric)),
    all(sapply(alpha, assertive::is_numeric)),
    identical(sort(names(abs_deviations)), sort(names(p))),
    identical(sort(names(abs_deviations)), sort(names(alpha))),
    msg = paste0("`abs_deviations`, `p`, and `alpha` must all be lists of ",
                 "numerics with the same named list elements")
  )

  assertthat::assert_that(
    all(names(abs_deviations) %in% possible_components),
    msg = paste0("named list elements must be part of '",
                 paste(possible_components, collapse = "', '"), "'.")
  )

  if ('survival' %in% names(abs_deviations)) {
    assertthat::assert_that(
      assertive::is_numeric(survival),
      all(data.table::between(survival, 0, 1)),
      msg = "`survival`, must be a numeric vector with all values between 0 and 1"
    )
  }

  # Calculate beta for all parameters other than survival -------------------

  basic_calculation <- function(parameter) {
    list(alpha = alpha[[parameter]],
         beta = calculate_beta(
           abs_deviations[[parameter]],
           p[[parameter]],
           alpha[[parameter]],
           log_transform = parameter %in% log_transform_components)
    )
  }

  parameters <- names(abs_deviations)
  hyperparameters <- lapply(setdiff(parameters, "survival"), function(parameter) {
    basic_calculation(parameter)
  })
  names(hyperparameters) <- setdiff(parameters, "survival")

  # Calculate beta for survival ---------------------------------------------

  # calculate beta for survival parameter
  calculate_death_lower_quantile <- function(alpha, beta, p, max_death_probability) {
    lower_quantile <- demUtils::logit(max_death_probability) -
      stats::qt(p = (1 + p) / 2, df = 2 * alpha) * sqrt(beta / alpha)

    lower_quantile <- demUtils::invlogit(lower_quantile)
    return(lower_quantile)
  }
  calculate_death_upper_quantile <- function(alpha, beta, p, max_death_probability) {
    upper_quantile <- demUtils::logit(max_death_probability) +
      stats::qt(p = (1 + p) / 2, df = 2 * alpha) * sqrt(beta / alpha)

    upper_quantile <- demUtils::invlogit(upper_quantile)
    return(upper_quantile)
  }
  survival_target <- function(absolute_deviation, max_death_probability) {
    beta <- calculate_beta(
      absolute_deviation,
      p[["survival"]],
      alpha[["survival"]],
      log_transform = F
    )
    lower_quantile <- calculate_death_lower_quantile(
      alpha[["survival"]],
      beta,
      p[["survival"]],
      max_death_probability
    )
    upper_quantile <- calculate_death_upper_quantile(
      alpha[["survival"]],
      beta,
      p[["survival"]],
      max_death_probability
    )

    if (upper_quantile <= lower_quantile) return(Inf)

    # TODO: why arbitrary cutoff?
    if (max_death_probability <= 0.5) {
      x <- (max_death_probability - lower_quantile) / max_death_probability
    } else {
      x <- (upper_quantile - max_death_probability) / max_death_probability
    }
    return((x - abs_deviations[["survival"]]) ^ 2)
  }

  if ("survival" %in% parameters) {
    min_surivival_proportion <- min(survival)
    abs_deviations_death_proportion <- stats::optimize(
      survival_target,
      interval = c(log(1 + abs_deviations[["survival"]]), log(2)),
      max_death_probability = 1 - min_surivival_proportion,
      tol = 1E-7
    )$minimum

    hyperparameters[["survival"]] <- list(
      alpha = alpha[["survival"]],
      beta = calculate_beta(
        abs_deviations_death_proportion,
        p[["survival"]],
        alpha[["survival"]],
        log_transform = F
      )
    )
  }

  return(hyperparameters)
}

#' @title Calculate the beta hyperparameter for popReconstruct variance
#'   parameters.
#'
#' @description Helper function for calculating the beta parameter of an inverse
#'   gamma distribution given absolute deviation from the popReconstruct
#'   elicitation statement. Used for quantifying measurement error in
#'   popReconstruct components.
#'
#' @param abs_deviation \[`numeric(1)`\]\cr
#'   Expert opinion about the absolute deviation \eqn{\eta_{\nu}} of the initial
#'   point estimates for a given probability interval \eqn{p_{\nu}}.
#' @param p \[`numeric(1)`\]\cr
#'   Central \eqn{p_{\nu}}% probability interval used in elicitation statement.
#' @param alpha \[`numeric(1)`\]\cr
#'   Alpha (shape) parameter of inverse gamma distribution.
#' @param log_transform \[`logical(1)`\]\cr
#'   Whether model component is log transformed in the popReconstruct model.
#'
#' @return \[`numeric(1)`\] Beta (scale) parameter of inverse gamma distribution
#' for the variance of the popReconstruct model component.
calculate_beta <- function(abs_deviation,
                           p,
                           alpha,
                           log_transform) {

  lower_quantile <- stats::qt(p = ((p + 1) / 2),
                              df = 2 * alpha)

  if (log_transform) {
    # 1 - abs_deviation is the more conservative choice
    beta <- alpha * ((log(1 - abs_deviation) / lower_quantile) ^ 2)
  } else {
    beta <- alpha * ((abs_deviation / lower_quantile) ^ 2)
  }
  return(beta)
}
