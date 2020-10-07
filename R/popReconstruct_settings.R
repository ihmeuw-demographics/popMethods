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
    if (identical(settings$sexes, "female")) settings$fixed_parameters <- "srb"
  }

  # add model components that are not being estimated
  model_components <- c("srb", "asfr", "baseline", "survival",
                        "net_migration", "immigration", "emigration")
  settings$fixed_parameters <- unique(c(settings$fixed_parameters,
                                        setdiff(model_components, names(inputs))))

  settings$estimated_parameters <- names(inputs)[!names(inputs) %in%
                                                   settings$fixed_parameters]
  settings$migration_parameters <- grep("migration$", names(inputs), value = T)

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
      ages = settings$ages_survival,
      ages_knots = settings$k_ages_survival,
      B_a = splines::bs(
        x = settings$ages_survival,
        knots = settings$k_ages_survival[-length(settings$k_ages_survival)],
        degree = 1
      ),
      transformation_name = "logit",
      transformation = demUtils::logit,
      inverse_transformation = demUtils::invlogit,
      fixed = "survival" %in% settings$fixed_parameters
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
      fixed = "emigration" %in% settings$fixed_parameters
    )
  )
  return(component_settings)
}
