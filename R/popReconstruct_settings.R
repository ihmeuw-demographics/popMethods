#' @title Validate or generate the optional settings used in the popReconstruct
#'   model
#'
#' @description
#' Adds general settings for the interval size 'int', each year
#' for which projected population estimates will exist 'years_projections', and
#' the years in which census data exists 'census_years.
#'
#' Also adds/validates settings related to specifying the B-spline linear basis
#' functions. For all the year-specific input components (all except for
#' 'baseline'), adds the 'k_years_{input_name}' setting which defines the start
#' of each calendar year interval that will be knots in the B-spline linear
#' basis functions. For all the age-specific input components (all except for
#' 'srb'), adds the 'k_ages_{input_name}' setting which defines the start of
#' each age-group interval that will be knots in the B-spline linear basis
#' functions.
#'
#' @inheritParams popReconstruct_fit
#'
#' @return Invisibly returns reference to modified `settings`.
create_optional_settings <- function(settings, inputs, data = NULL) {

  # age and year intervals
  settings$int <- unique(diff(settings$ages))

  settings$years_projections <- c(settings$years,
                                  max(settings$years) + settings$int)
  if (!is.null(data)) {
    settings$years_census <- unique(data$population$year)
  }

  # "fixed_parameters" setting
  if ("fixed_parameters" %in% names(settings)) {
    assertthat::assert_that(
      all(settings$fixed_parameters %in% names(inputs)),
      msg = "setting 'fixed_parameters' must correspond to given `inputs`"
    )
  } else {
    settings$fixed_parameters <- NULL
  }
  # if not estimating both sexes don't estimate srb
  if (identical(settings$sexes, "female")) {
    settings$fixed_parameters <- c(settings$fixed_parameters, "srb")
  }

  # add model components that are not being estimated
  mortality_parameters <- c("survival", "mx", "non_terminal_ax", "terminal_ax")
  migration_parameters <- c("net_migration", "immigration", "emigration")
  model_components <- c(
    "srb", "asfr", "baseline",
    mortality_parameters, migration_parameters
  )
  settings$fixed_parameters <- unique(c(settings$fixed_parameters,
                                        setdiff(model_components, names(inputs))))
  settings$estimated_parameters <- names(inputs)[!names(inputs) %in%
                                                   settings$fixed_parameters]
  settings$mortality_parameters <- mortality_parameters[mortality_parameters %in% names(inputs)]
  settings$migration_parameters <- migration_parameters[migration_parameters %in% names(inputs)]

  # "n_draws" setting
  if ("n_draws" %in% names(settings)) {
    assertthat::assert_that(
      assertthat::is.number(settings$n_draws),
      settings$n_draws > 0,
      msg = "setting 'n_draws' must be a number greater than zero"
    )
  } else {
    settings$n_draws <- 1000
  }

  all_components <- c(settings$estimated_parameters, settings$fixed_parameters)

  # "k_years_{input_name}" setting
  for (comp in setdiff(all_components, "baseline")) {
    setting_name <- paste0("k_years_", comp)
    if (setting_name %in% names(settings)) {
      setting_value <- settings[[setting_name]]
      assertthat::assert_that(
        assertive::is_numeric(setting_value),
        all(setting_value %in% settings$years),
        msg = paste0("setting '", setting_name, "' must be a numeric ",
                     "corresponding to the start of calendar year intervals ",
                     "as specified in the setting 'years'")
      )
    } else {
      settings[[setting_name]] <- settings$years
    }
  }

  # "k_ages_{input_name}" setting
  for (comp in setdiff(all_components, "srb")) {
    setting_name <- paste0("k_ages_", comp)

    # determine the expected ages for each component
    expected_ages <- settings$ages
    if (paste0("ages_", comp) %in% names(settings)) {
      expected_ages <- settings[[paste0("ages_", comp)]]
    }

    if (setting_name %in% names(settings)) {
      setting_value <- settings[[setting_name]]

      assertthat::assert_that(
        assertive::is_numeric(setting_value),
        all(setting_value %in% expected_ages),
        msg = paste0("setting '", setting_name, "' must be a numeric ",
                     "corresponding to the start of age-group intervals ",
                     "as specified in the setting 'ages' or ",
                     "'ages_{input_name}'")
      )
    } else {
      settings[[setting_name]] <- expected_ages
    }
  }

  return(invisible(settings))
}

#' @title Helper function to create a detailed list of settings specific to
#'   popReconstruct ccmpp inputs and outputs.
#'
#' @inheritParams popReconstruct_fit
#'
#' @section Detailed settings:
#' * measure_type: whether the component is either stock (measured at a specific
#'     point in time) or flow (measured over a period of time). This determines
#'     if the component will have a 'year' or 'year_start' and 'year_end'
#'     columns.
#' * id_cols: the id columns defining the input and output \[`data.table()`\].
#' * years: the 'year_start' or 'year' values the input/output will have.
#' * years_knots: the start of each calendar year interval that will be knots
#'     in the B-spline linear basis functions.
#' * B_t: the B-spline basis matrix for a linear spline with knots specified at
#'     'years_knots'.
#' * sexes: the 'sex' values the input/output will have.
#' * ages: the 'age_start' values the input/output will have.
#' * age_knots: the start of each age group interval that will be knots in the
#'     B-spline linear basis functions.
#' * B_a: the B-spline basis matrix for a linear spline with knots specified at
#'     'ages_knots'.
#' * transformation_name: character name of the transformation.
#' * transformation: the transformation used when modeling each component.
#' * inverse_transformation: the inverse transformation used when modeling each
#'     component.
#' * transformation_arguments: additional arguments to pass to the 'transformation'
#'     and 'inverse_transformation' functions.
#' * fixed: whether the ccmpp input component is to be estimated or fixed in the
#'     model.
#'
#' @return \[`list()`\] of detailed settings for each ccmpp input and
#' 'population'.
create_detailed_settings <- function(settings) {

  component_settings <- list(
    population = list(
      measure_type = "stock",
      id_cols = c("year", "sex", "age_start", "age_end"),
      years = settings$years_projections,
      sexes = settings$sexes,
      ages = settings$ages,
      transformation_name = "log",
      transformation = log,
      inverse_transformation = exp,
      transformation_arguments = NULL,
      fixed = FALSE
    ),
    srb = list(
      measure_type = "flow",
      id_cols = c("year_start", "year_end"),
      years = settings$years,
      years_knots = settings$k_years_srb,
      B_t = splines::bs(
        x = settings$years,
        knots = settings$k_years_srb[-length(settings$k_years_srb)],
        degree = 1
      ),
      transformation_name = "log",
      transformation = log,
      inverse_transformation = exp,
      transformation_arguments = NULL,
      fixed = "srb" %in% settings$fixed_parameters
    ),
    asfr = list(
      measure_type = "flow",
      id_cols = c("year_start", "year_end", "age_start", "age_end"),
      years = settings$years,
      years_knots = settings$k_years_asfr,
      B_t = splines::bs(
        x = settings$years,
        knots = settings$k_years_asfr[-length(settings$k_years_asfr)],
        degree = 1
      ),
      ages = settings$ages_asfr,
      ages_knots = settings$k_ages_asfr,
      B_a = splines::bs(
        x = settings$ages_asfr,
        knots = settings$k_ages_asfr[-length(settings$k_ages_asfr)],
        degree = 1
      ),
      transformation_name = "log",
      transformation = log,
      inverse_transformation = exp,
      transformation_arguments = NULL,
      fixed = "asfr" %in% settings$fixed_parameters
    ),
    baseline = list(
      measure_type = "stock",
      id_cols = c("year", "sex", "age_start", "age_end"),
      years = min(settings$years),
      sexes = settings$sexes,
      ages = settings$ages,
      ages_knots = settings$k_ages_baseline,
      B_a = splines::bs(
        x = settings$ages,
        knots = settings$k_ages_baseline[-length(settings$k_ages_baseline)],
        degree = 1
      ),
      transformation_name = "log",
      transformation = log,
      inverse_transformation = exp,
      transformation_arguments = NULL,
      fixed = "baseline" %in% settings$fixed_parameters
    ),
    survival = list(
      measure_type = "flow",
      id_cols = c("year_start", "year_end", "sex", "age_start", "age_end"),
      years = settings$years,
      years_knots = settings$k_years_survival,
      B_t = splines::bs(
        x = settings$years,
        knots = settings$k_years_survival[-length(settings$k_years_survival)],
        degree = 1
      ),
      sexes = settings$sexes,
      ages = settings$ages_mortality,
      ages_knots = settings$k_ages_survival,
      B_a = splines::bs(
        x = settings$ages_mortality,
        knots = settings$k_ages_survival[-length(settings$k_ages_survival)],
        degree = 1
      ),
      transformation_name = "logit",
      transformation = demUtils::logit,
      inverse_transformation = demUtils::invlogit,
      transformation_arguments = NULL,
      fixed = "survival" %in% settings$fixed_parameters
    ),
    mx = list(
      measure_type = "flow",
      id_cols = c("year_start", "year_end", "sex", "age_start", "age_end"),
      years = settings$years,
      years_knots = settings$k_years_mx,
      B_t = splines::bs(
        x = settings$years,
        knots = settings$k_years_mx[-length(settings$k_years_mx)],
        degree = 1
      ),
      sexes = settings$sexes,
      ages = settings$ages_mortality,
      ages_knots = settings$k_ages_mx,
      B_a = splines::bs(
        x = settings$ages_mortality,
        knots = settings$k_ages_mx[-length(settings$k_ages_mx)],
        degree = 1
      ),
      transformation_name = "log",
      transformation = log,
      inverse_transformation = exp,
      transformation_arguments = NULL,
      fixed = "mx" %in% settings$fixed_parameters
    ),
    non_terminal_ax = list(
      measure_type = "flow",
      id_cols = c("year_start", "year_end", "sex", "age_start", "age_end"),
      years = settings$years,
      years_knots = settings$k_years_non_terminal_ax,
      B_t = splines::bs(
        x = settings$years,
        knots = settings$k_years_non_terminal_ax[-length(settings$k_years_non_terminal_ax)],
        degree = 1
      ),
      sexes = settings$sexes,
      ages = settings$ages_mortality[-length(settings$ages_mortality)],
      ages_knots = settings$k_ages_non_terminal_ax,
      B_a = splines::bs(
        x = settings$ages,
        knots = settings$k_ages_non_terminal_ax[-length(settings$k_ages_non_terminal_ax)],
        degree = 1
      ),
      # ax (average years lived by those dying in the age interval) for all age
      # groups except the terminal age group must be constrained to be between 0
      # and the age interval.
      transformation_name = "bounded_logit",
      transformation = demUtils::logit,
      inverse_transformation = demUtils::invlogit,
      transformation_arguments = list(domain_lower = 0, domain_upper = settings$int),
      fixed = "non_terminal_ax" %in% settings$fixed_parameters
    ),
    terminal_ax = list(
      measure_type = "flow",
      id_cols = c("year_start", "year_end", "sex", "age_start", "age_end"),
      years = settings$years,
      years_knots = settings$k_years_terminal_ax,
      B_t = splines::bs(
        x = settings$years,
        knots = settings$k_years_terminal_ax[-length(settings$k_years_terminal_ax)],
        degree = 1
      ),
      sexes = settings$sexes,
      ages = settings$ages_mortality[length(settings$ages_mortality)],
      ages_knots = settings$ages_mortality[length(settings$ages_mortality)],
      B_a = NULL,
      transformation_name = "log",
      transformation = log,
      inverse_transformation = exp,
      transformation_arguments = NULL,
      fixed = "terminal_ax" %in% settings$fixed_parameters
    ),
    net_migration = list(
      measure_type = "flow",
      id_cols = c("year_start", "year_end", "sex", "age_start", "age_end"),
      years = settings$years,
      years_knots = settings$k_years_net_migration,
      B_t = splines::bs(
        x = settings$years,
        knots = settings$k_years_net_migration[-length(settings$k_years_net_migration)],
        degree = 1
      ),
      sexes = settings$sexes,
      ages = settings$ages,
      ages_knots = settings$k_ages_net_migration,
      B_a = splines::bs(
        x = settings$ages,
        knots = settings$k_ages_net_migration[-length(settings$k_ages_net_migration)],
        degree = 1
      ),
      transformation_name = NULL,
      transformation = NULL,
      inverse_transformation = NULL,
      transformation_arguments = NULL,
      fixed = "net_migration" %in% settings$fixed_parameters
    ),
    immigration = list(
      measure_type = "flow",
      id_cols = c("year_start", "year_end", "sex", "age_start", "age_end"),
      years = settings$years,
      years_knots = settings$k_years_immigration,
      B_t = splines::bs(
        x = settings$years,
        knots = settings$k_years_immigration[-length(settings$k_years_immigration)],
        degree = 1
      ),
      sexes = settings$sexes,
      ages = settings$ages,
      ages_knots = settings$k_ages_immigration,
      B_a = splines::bs(
        x = settings$ages,
        knots = settings$k_ages_immigration[-length(settings$k_ages_immigration)],
        degree = 1
      ),
      transformation_name = "log",
      transformation = log,
      inverse_transformation = exp,
      transformation_arguments = NULL,
      fixed = "immigration" %in% settings$fixed_parameters
    ),
    emigration = list(
      measure_type = "flow",
      id_cols = c("year_start", "year_end", "sex", "age_start", "age_end"),
      years = settings$years,
      years_knots = settings$k_years_emigration,
      B_t = splines::bs(
        x = settings$years,
        knots = settings$k_years_emigration[-length(settings$k_years_emigration)],
        degree = 1
      ),
      sexes = settings$sexes,
      ages = settings$ages,
      ages_knots = settings$k_ages_emigration,
      B_a = splines::bs(
        x = settings$ages,
        knots = settings$k_ages_emigration[-length(settings$k_ages_emigration)],
        degree = 1
      ),
      transformation_name = "log",
      transformation = log,
      inverse_transformation = exp,
      transformation_arguments = NULL,
      fixed = "emigration" %in% settings$fixed_parameters
    )
  )
  return(component_settings)
}

#' @title Helper function to transform values in a data.table
#'
#' @param dt \[`data.table(1)`\]\cr
#'   The data.table to apply the transformation functions to.
#' @param value_col \[`character(1)`\]\cr
#'   The column in `dt` to transform.
#' @param transformation \[`function(1)`\]\cr
#'   NULL if no transformation. Otherwise transformation function to
#'   apply to `value_col`.
#' @param transformation_arguments \[`list()`\]\cr
#'   NULL if no transformation. Otherwise list of arguments to provide to the
#'   corresponding `transformation` function(s).
#'
#' @return invisibly return `dt` with transformed values.
#'
#' @examples
#' popMethods:::transform_dt(
#'   dt = data.table::data.table(year = 1950, value = 2.5),
#'   value_col = "value",
#'   transformation = demUtils::logit,
#'   transformation_arguments = list(domain_lower = 0, domain_upper = 5)
#' )
transform_dt <- function(dt,
                         value_col,
                         transformation,
                         transformation_arguments) {

  if (!is.null(transformation)) {
    setnames(dt, value_col, "initial_value")
    dt[
      ,
      transformed_value := do.call(
        what = transformation,
        args = c(list(x = initial_value), transformation_arguments)
      )]
    setnames(dt, "transformed_value", value_col)
    dt[, initial_value := NULL]
  }
  return(dt)
}
