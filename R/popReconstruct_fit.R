#' @title Implement the Bayesian hierarchical PopReconstruct model
#'
#' @description
#' The popReconstruct model embeds the cohort component method of population
#' projection (ccmpp) in a Bayesian hierarchical model. It reconciles
#' ccmpp input preliminary estimates of sex-ratio at birth, fertility,
#' mortality, migration, and a baseline population with population data while
#' accounting for measurement error in these initial estimates.
#'
#' `popReconstruct_fit` fits the model using [Stan](https://mc-stan.org/) or
#' [TMB](https://kaskr.github.io/adcomp/Introduction.html).
#'
#' `popReconstruct_posterior_draws` produces draws from the posterior
#' distribution.
#'
#' `popReconstruct_prior_draws` produces draws from the prior distribution.
#'
#' `popReconstruct_summarize_draws` produces summary statistics of draws from
#' the popReconstruct model using `demUtils::summarize_dt`. The `summarize_cols`
#' parameter should include 'chain', 'chain_draw' and 'draw' when the model was
#' fit with stan and just 'draw' when the model was fit with tmb.
#'
#' @inheritParams demCore::ccmpp
#' @param data \[`list()`\]\cr
#'   \[`data.table()`\] for each type of data input for the model (currently
#'   only 'population'). See **Section: Data** for more information.
#' @param hyperparameters \[`list()`\]\cr
#'   List for each model component ('srb', 'asfr', ...) that itself contains a
#'   list specifying the alpha (shape) and beta (scale) hyperparameters for the
#'   variance of each model component.
#' @param settings \[`list()`\]\cr
#'   Named list of settings for running the popReconstruct model with. The
#'   required settings are the same as those required for `demCore::ccmpp()`,
#'   see **Section: Settings** for these. The optional settings specific to
#'   popReconstruct are described in
#'   **Section: Optional popReconstruct Settings**
#' @param software \[`character()`\]\cr
#'   Statistical modeling software to fit popReconstruct with. Either 'stan' or
#'   'tmb' are available.
#'
#' @param fit \[`stanfit(1)`\] or \[`sdreport(1)`\]\cr
#'   Model fit object returned by `popReconstruct_fit`.
#' @param method_name \[`character(1)`\]\cr
#'   Description to assign to column 'method' in each returned data.table.
#' @param chunk_size \[`integer(1)`\]\cr
#'   number of draws to sample from the prior at once.
#'
#' @param draws \[`list()`\] of \[`data.table()`\]\cr
#'   Draws from the posterior distribution of the popReconstruct model as
#'   returned by `popReconstruct_posterior_draws`.
#' @inheritParams demUtils::summarize_dt
#'
#' @param ... For `popReconstruct_fit` additional arguments to pass to
#'   `rstan::sampling` if fitting the model with 'stan' or `optimx::optimx` if
#'   fitting the model with 'tmb'. For `popReconstruct_summarize_draws`
#'   additional arguments to pass to `demUtils::summarize_dt`.
#'
#' @inheritSection demCore::ccmpp Settings
#'
#' @section Optional popReconstruct Settings:
#'   * fixed_parameters: \[`character()`\]\cr
#'     Names of `ccmpp` input components that should be fixed in the model,
#'     corresponds to the names of the `input` list. Defaults to empty if
#'     fitting the model for both sexes, if only 'female' then defaults to
#'     'srb'.
#'   * n_draws: \[`character()`\]\cr
#'     The number of draws to take from the output posterior distribution. If
#'     fitting the model with stan then the number of draws is dictated by
#'     `chains`, `iterations`, `burn-in`, `thin` arguments to `rstan::sampling`;
#'     `n_draws` can then be used to subsample down.
#'   * k_years_{input_name}: \[`integer()`\]\cr
#'     The start of each calendar year interval that will be knots in the
#'     B-spline linear basis functions for the year-specific input estimates,
#'     so not relevant for 'baseline'. Defaults to the setting for 'years' which
#'     is equal to the original popReconstruct model without linear B-splines
#'     to reduce year-specific parameters.
#'   * k_ages_{input_name}: \[`integer()`\]\cr
#'     The start of each age-group interval that will be knots in the B-spline
#'     linear basis functions for the age-specific input estimates, so not
#'     relevant for 'srb'. Defaults to the setting for 'ages' (or
#'     'ages_{input_name}') which is equal to the original popReconstruct model
#'     without linear B-splines to reduce age-specific parameters.
#'
#' @inheritSection demCore::ccmpp Inputs
#'
#' @section `popReconstruct_fit` Value:
#' If fit using 'stan' an object of class `stanfit` and if fit using 'tmb' an
#' object of class `sdreport`. Either represents the fitted results that can be
#' extracted using `popReconstruct_posterior_draws()`. Stan has helper packages
#' that can be used to explore the model fit through `shinystan`, `bayesplot`,
#' etc.
#'
#' @section `popReconstruct_posterior_draws` and `popReconstruct_prior_draws` Value:
#'
#' `popReconstruct_posterior_draws` and `popReconstruct_prior_draws` return a
#' named \[`list()`\] of \[`data.table()`\] for draws from the posterior and
#' prior distribution respectively for each ccmpp input component along with the
#' associated offset and spline offset parameters. Draws for the 'variance' and
#' projected 'population' are also included.
#'
#' The returned \[`data.table()`\]s will include columns related to the draw
#' number:
#' * draw: the draw index.
#' * chain: if using stan to fit the model, the chain the draw came from.
#' * chain_draw: if using stan to fit the model, the draw index for that chain.
#'
#' population: \[`data.table()`\]\cr
#'   * year: \[`integer()`\] mid-year for population estimate.
#'   Corresponds to 'years' `setting` plus an extra interval.
#'   * sex: \[`character()`\] either 'female' or 'male'. Corresponds to 'sexes'
#'   `setting`.
#'   * age_start: \[`integer()`\] start of the age group (inclusive).
#'   Corresponds to 'ages' `setting`.
#'   * age_end: \[`integer()`\] end of the age group (exclusive).
#'   * draw: \[`integer()`\] draw index.
#'   * `value_col`: \[`numeric()`\] projected population count estimates.
#'
#' variance: \[`data.table()`\]\cr
#'   * parameter: \[`character()`\] the ccmpp component the variance term
#'   corresponds to.
#'   * draw: \[`integer()`\] draw index.
#'   * `value_col`: \[`numeric()`\] estimated variance value..
#'
#' For each ccmpp component ('srb', 'asfr', etc.) there will be three
#' \[`data.table()`\] outputs included in the returned \[`list()`\]. These will
#' each have the same columns as described in the **Section: Inputs** but with
#' additional draw related columns as described above.
#'
#' 1. The offset parameters representing the values of the piecewise linear
#' function at the year- and age-specific knots. These are the estimated
#' deviations from the initial ccmpp input estimates at the exact knots. The
#' estimated values are in the transformed modeling space (log transformed for
#' 'srb', logit transformed for 'survival', etc.).
#' 2. The spline offset values calculated by multiplying the fixed B-spline
#' linear basis functions by the offset parameters. These are the estimated
#' deviations from the initial ccmpp input estimates at all years and ages. The
#' estimated values are in the transformed modeling space.
#' 3. The actual ccmpp input posterior draws after combining the spline offset
#' values with the initial ccmpp input estimates and after applying the inverse
#' of the transformation used to model each component.
#'
#' @section `popReconstruct_summarize_draws` Value:
#' Returns a named \[`list()`\] of \[`data.table()`\] for summary statistics
#' of the input `draws`. Each [data.table()] will have the `id_cols` for each
#' component (minus the `summarize_cols`) plus summary statistic columns.
#' The summary statistic columns have the same name as each function specified
#' in `summary_fun` and the quantiles are named like 'q_(`probs` * 100)'. See
#' `demUtils::summarize_dt` for more information.
#'
#' @references
#' Wheldon, Mark C., Adrian E. Raftery, Samuel J. Clark, and Patrick Gerland.
#' 2013. “Reconstructing Past Populations With Uncertainty From Fragmentary
#' Data.” Journal of the American Statistical Association 108 (501): 96–110.
#' [https://doi.org/10.1080/01621459.2012.737729](https://doi.org/10.1080/01621459.2012.737729).
#'
#' [popReconstruct R Package](https://cran.r-project.org/web/packages/popReconstruct/popReconstruct.pdf)
#'
#' Wheldon, Mark C., Adrian E. Raftery, Samuel J. Clark, and Patrick Gerland.
#' 2015. “Bayesian Reconstruction of Two-Sex Populations by Age: Estimating Sex
#' Ratios at Birth and Sex Ratios of Mortality.” Journal of the Royal
#' Statistical Society. Series A: Statistics in Society 178 (4): 977–1007.
#' [https://doi.org/10.1111/rssa.12104](https://doi.org/10.1111/rssa.12104).
#'
#' [markalava/Bayesian-Reconstruction github repo](https://github.com/markalava/Bayesian-Reconstruction)
#'
#' @examples
#' # specify settings for this example
#' settings = list(
#'   years = seq(1960, 2000, 5),
#'   sexes = c("female"),
#'   ages = seq(0, 80, 5),
#'   ages_survival = seq(0, 85, 5),
#'   ages_asfr = seq(15, 45, 5),
#'
#'   n_draws = 10
#' )
#'
#' # hyperparameters for the variance prior distribution, represents measurement error
#' hyperparameters <- list(asfr = list(alpha = 1, beta = 0.0109),
#'                         population = list(alpha = 1, beta = 0.0109),
#'                         survival = list(alpha = 1, beta = 0.0109),
#'                         net_migration = list(alpha = 1, beta = 0.0436))
#'
#' fit_stan <- popMethods::popReconstruct_fit(
#'   inputs = demCore::burkina_faso_initial_estimates,
#'   data = demCore::burkina_faso_data,
#'   hyperparameters = hyperparameters,
#'   settings = settings,
#'   value_col = "value",
#'   software = "stan",
#'   chains = 1,
#'   iter = 200,
#'   warmup = 100,
#'   thin = 2
#' )
#'
#' draws_stan <- popMethods::popReconstruct_posterior_draws(
#'   fit = fit_stan,
#'   inputs = demCore::burkina_faso_initial_estimates,
#'   settings = settings,
#'   value_col = "value",
#'   software = "stan",
#'   method_name = "original"
#' )
#'
#' summary_stan <- popMethods::popReconstruct_summarize_draws(
#'   draws = draws_stan
#' )
#'
#' @seealso demCore::ccmpp
#' @seealso demUtils::summarize_dt
#' @family popReconstruct
#'
#' @import Rcpp
#' @import methods
#' @rawNamespace useDynLib(popMethods, .registration = TRUE)
#' @rawNamespace useDynLib(popMethods_TMBExports)
#'
#' @export
#' @rdname popReconstruct
popReconstruct_fit <- function(inputs,
                               data,
                               hyperparameters,
                               settings,
                               value_col,
                               software,
                               ...) {

  # Validate arguments ------------------------------------------------------

  # check `software` argument
  assertthat::assert_that(
    assertthat::is.string(software),
    software %in% c("stan", "tmb"),
    msg = "`software` must be 'stan' or 'tmb'."
  )

  demCore:::validate_ccmpp_inputs(inputs, settings, value_col)

  # add optional settings needed for popReconstruct
  settings <- copy(settings)
  settings <- create_optional_settings(settings, inputs, data)
  detailed_settings <- create_detailed_settings(settings)

  validate_popReconstruct_hyperparameters(hyperparameters, inputs, settings)
  validate_popReconstruct_data(data, settings, value_col)

  # Create input objects for fitting model ----------------------------------

  # determine columns of projected population matrix in Stan that corresponds
  # to census years
  index <- 1:length(settings$years_projections)
  if (software == "tmb") index <- index - 1
  index_census_years <- index[settings$years_projections %in%
                                settings$years_census]

  # determine number of younger age groups that are not included in the
  # reproductive ages
  A_f_offset = sum(settings$ages < min(settings$ages_asfr))

  # prepare input data for the model
  input_parameters <- list()
  input_data <- list(
    sexes = length(settings$sexes),
    interval = settings$int,
    A = length(settings$ages),
    A_f = length(settings$ages_asfr),
    A_f_offset = A_f_offset,
    Y = length(settings$years),
    Y_p = length(settings$years_census),
    pop_data_years_index = index_census_years
  )

  prep_input_array <- function(component, list_dt, detailed_settings) {
    dt <- copy(list_dt[[component]])

    # transform value
    setnames(dt, value_col, "value")
    transformation <- detailed_settings[[component]][["transformation"]]
    if (!is.null(transformation)) {
      dt <- dt[, value := transformation(value)]
    }

    # convert from dt to matrix format
    mdt <- demCore:::dt_to_matrix(
      dt = dt,
      id_cols = setdiff(names(dt), "value"),
      value_col = "value"
    )

    # if sex-specific format as array
    adt <- mdt
    if (assertive::is_list(mdt)) {
      n_sexes <- length(mdt)
      n_ages <- nrow(mdt[[1]])
      n_years <- ncol(mdt[[1]])
      if (software == "stan") {
        adt <- array(dim = c(n_sexes, n_ages, n_years))
        for (s in 1:n_sexes) {
          adt[s, , ] <- mdt[[s]]
        }
      } else if (software == "tmb") {
        adt <- array(
          unlist(mdt),
          dim = c(n_ages, n_years, n_sexes)
        )
      }
    }

    return(adt)
  }

  # prepare data inputs
  for (component in c("population")) {

    # add initial input estimates
    data_adt <- prep_input_array(component, data, detailed_settings)
    transformation_name <- detailed_settings[[component]][["transformation_name"]]
    input_name <- paste0("input_", transformation_name,
                         if (!is.null(transformation_name)) "_",
                         component, "_data")
    input_data[[input_name]] <- data_adt
  }

  # prepare initial estimates and parameters for all possible components
  # specified in the stan and tmb models
  all_ccmpp_inputs <- c("srb", "asfr", "baseline", "survival", "net_migration")
  for (component in all_ccmpp_inputs) {

    sexes <- detailed_settings[[component]][["sexes"]]
    years_knots <- detailed_settings[[component]][["years_knots"]]
    ages_knots <- detailed_settings[[component]][["ages_knots"]]

    # add B-spline basis functions
    if ("B_t" %in% names(detailed_settings[[component]])) {
      input_data[[paste0("N_k_t_", component)]] <- length(years_knots)
      input_data[[paste0("B_t_", component)]] <-
        detailed_settings[[component]][["B_t"]]
    }
    if ("B_a" %in% names(detailed_settings[[component]])) {
      input_data[[paste0("N_k_a_", component)]] <- length(ages_knots)
      input_data[[paste0("B_a_", component)]] <-
        detailed_settings[[component]][["B_a"]]
    }

    # add initial input estimates
    input_adt <- prep_input_array(component, inputs, detailed_settings)
    transformation_name <- detailed_settings[[component]][["transformation_name"]]
    input_name <- paste0("input_", transformation_name,
                         if (!is.null(transformation_name)) "_",
                         component)
    input_data[[input_name]] <- input_adt

    # add initial offset parameter values for the model
    offset_dimension <- c(
      ifelse(!is.null(ages_knots), length(ages_knots), 1),
      ifelse(!is.null(years_knots), length(years_knots), 1)
    )
    if (software == "stan") {
      # add sex dimension to beginning for stan
      if (!is.null(sexes)) {
        offset_dimension <- c(length(sexes), offset_dimension)
      }

      # add indicator variable for whether component is being estimated or not
      fixed <- detailed_settings[[component]][["fixed"]]
      input_data[[paste0("estimate_", component)]] <- as.integer(!fixed)
      offset_dimension <- c(input_data[[paste0("estimate_", component)]],
                            offset_dimension)

    } else if (software == "tmb") {
      # add sex dimension to end for tmb
      if (!is.null(sexes)) {
        offset_dimension <- c(offset_dimension, length(sexes))
      }
    }
    offset_name <- paste0("offset_", transformation_name,
                          if (!is.null(transformation_name)) "_",
                          component)
    input_parameters[[offset_name]] <- array(0, dim = offset_dimension)
  }

  # prepare data and parameters for all possible components specified in the
  # stan and tmb models
  all_variance_inputs <- c("srb", "asfr", "population", "survival",
                           "net_migration")
  for (component in all_variance_inputs) {
    # add alpha and beta hyperparameters
    for (param in c("alpha", "beta")) {
      value <- 1
      if (component %in% names(hyperparameters)) {
        value <- hyperparameters[[component]][[param]]
      }
      input_data[[paste0(param, "_", component)]] <- value
    }

    # tmb requires an initial value for all parameters
    if (software == "tmb") {
      initial_variance <- 0.02 # TODO: make an option?
      variance_name <- paste0("log_sigma2_", component)
      input_parameters[[variance_name]] <- log(initial_variance)
    }
  }

  # actually fit the popReconstruct model using the specified software
  if (software == "stan") {

    # stan requires an initial parameter value for each chain
    # parameters without specified initial values will be randomly chosen
    chains <- ifelse("chains" %in% names(list(...)), list(...)[["chains"]], 2)
    input_parameters <- lapply(1:chains, function(i) input_parameters)

    # draw samples from the model
    fit <- rstan::sampling(
      object = stanmodels$popReconstruct,
      data = input_data,
      init = input_parameters,
      ...
    )

  } else if (software == "tmb") {

    # check if fixing certain parameters
    map <- NULL
    if (length(settings$fixed_parameters) > 0) {
      fix_params <- input_parameters[grepl(settings$fixed_parameters,
                                           names(input_parameters))]
      map <- lapply(names(fix_params), function(p) {
        component_map <- fix_params[[p]]
        if (assertthat::is.number(fix_params[[p]])) { # log_sigma2 parameters
          component_map <- factor(NA)
        } else { # ccmpp input parameters
          component_map <-  rep(factor(NA), length(fix_params[[p]]))
          dim(component_map) <- dim(fix_params[[p]])
        }
        return(component_map)
      })
      names(map) <- names(fix_params)
    }

    # make objective function
    input_data <- c(list(model = "popReconstruct"), input_data)
    obj <- TMB::MakeADFun(
      data = input_data,
      input_parameters,
      DLL = "popMethods_TMBExports",
      random = grep("^offset_", names(input_parameters), value = T),
      map = map
    )

    # optimize objective function
    opt <- optimx::optimx(
      par = obj$par,
      fn = function(x) as.numeric(obj$fn(x)),
      gr = obj$gr,
      method = "nlminb",
      ...
    )
    if (opt$convcod != 0) stop("tmb popReconstruct model did not converge")

    # calculate standard deviations of all model parameters
    fit <- TMB::sdreport(obj, getJointPrecision = T)
  }
  return(fit)
}

validate_popReconstruct_hyperparameters <- function(hyperparameters,
                                                    inputs,
                                                    settings) {

  possible_components <- c("srb", "asfr", "population", "survival",
                           settings$migration_parameters)
  expected_components <- setdiff(possible_components, settings$fixed_parameters)

  assertthat::assert_that(
    class(hyperparameters) == "list",
    identical(sort(names(hyperparameters)), sort(expected_components)),
    all(sapply(hyperparameters, class) == "list"),
    all(sapply(hyperparameters, function(comp) sapply(comp, assertthat::is.number))),
    all(sapply(hyperparameters, function(comp) identical(sort(names(comp)),
                                                         c("alpha", "beta")))),
    msg = paste0("`hyperparameters` must be a list of lists with the first ",
                 "level containing each non-fixed model component and the ",
                 "second level containing named 'alpha' and 'beta' parameters.")
  )

}

validate_popReconstruct_data <- function(data, settings, value_col) {

  # check all required columns are present in `data`
  component_cols <- list(
    "population" = c("year", "sex", "age_start", "age_end")
  )

  component_ids <- list(
    "population" = list(year = settings$years_projections,
                             sex = settings$sexes,
                             age_start = settings$ages)
  )
  assertthat::assert_that(
    class(data) == "list",
    names(data) %in% names(component_ids),
    all(mapply(assertive::is_data.table, data)),
    msg = paste0("`data` must be a list of data.tables with named elements for '",
                 paste(names(component_ids), collapse = "', '"), "'")
  )

  for (component in names(component_ids)) {
    required_cols <- component_cols[[component]]
    assertthat::assert_that(
      all(required_cols %in% names(data[[component]])),
      value_col %in% names(data[[component]]),
      msg = paste0(component, " must be included in `data` with columns '",
                   paste(c(required_cols, value_col), collapse = "', '"), "'.")
    )

    ids <- component_ids[[component]]
    for (col in names(ids)) {
      assertthat::assert_that(
        all(data[[component]][[col]] %in% ids[[col]]),
        msg = paste0(component, " of `data` must contain '", col, "' values ",
                     "that are part of '", paste(ids[[col]], collapse = "', '"),"'")
      )
    }
    assertable::assert_values(data[[component]], colnames = value_col,
                              test = "not_na", quiet = T)
  }

}
