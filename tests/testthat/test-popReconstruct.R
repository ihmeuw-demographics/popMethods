
# load stan helper functions for testing
rstan::expose_stan_functions(stanmodel = stanmodels$popReconstruct)

set.seed(1)

settings = list(
  years = seq(1960, 2000, 5),
  sexes = c("female"),
  ages = seq(0, 80, 5),
  ages_mortality = seq(0, 85, 5),
  ages_asfr = seq(15, 45, 5),

  n_draws = 100
)
hyperparameters <- list(
  asfr = list(alpha = 1, beta = 0.0109),
  population = list(alpha = 1, beta = 0.0109),
  survival = list(alpha = 1, beta = 0.0109),
  net_migration = list(alpha = 1, beta = 0.0436)
)

# popReconstruct testing helper functions ---------------------------------

test_combinations <- function(desc,
                              inputs,
                              data,
                              hyperparameters,
                              settings,
                              software,
                              count_parameters,
                              test_software = c("stan", "tmb", "prior")) {

  if ("stan" %in% test_software) {
    description <- paste(desc, "(stan)")
    testthat::test_that(description, {
      test_fit(
        inputs = inputs,
        data = data,
        hyperparameters = hyperparameters,
        settings = settings,
        software = "stan",
        count_parameters = count_parameters,
        chains = 1, warmup = 100, iter = 200, thin = 2, seed = 3
      )
    })
  }

  if ("tmb" %in% test_software) {
    description <- paste(desc, "(tmb)")
    testthat::test_that(description, {
      test_fit(
        inputs = inputs,
        data = data,
        hyperparameters = hyperparameters,
        settings = settings,
        software = "tmb",
        count_parameters = count_parameters
      )
    })
  }

  if ("prior" %in% test_software) {
    description <- paste(desc, "(prior)")
    testthat::test_that(description, {
      test_prior(
        inputs = inputs,
        hyperparameters = hyperparameters,
        settings = settings,
        count_parameters = count_parameters
      )
    })
  }
}

test_fit <- function(inputs,
                     data,
                     hyperparameters,
                     settings,
                     software,
                     count_parameters,
                     ...) {

  fit <- testthat::expect_error(
    popMethods::popReconstruct_fit(
      inputs = inputs,
      data = data,
      hyperparameters = hyperparameters,
      settings = settings,
      value_col = "value",
      software = software,
      ...
    ),
    NA
  )

  draws <- testthat::expect_silent(
    popMethods::popReconstruct_posterior_draws(
      fit = fit,
      inputs = inputs,
      settings = settings,
      value_col = "value",
      software = software,
      method_name = "original"
    )
  )
  testthat::expect_equal("list", class(draws))
  testthat::expect_equal("data.table", unique(sapply(draws, class)[1, ]))

  testthat::expect_silent(
    draws <- popMethods::popReconstruct_count_space_parameters(
      draws = draws,
      settings = settings,
      parameters = count_parameters,
      value_col = "value",
      quiet = TRUE
    )
  )
  testthat::expect_true(all(count_parameters %in% names(draws)))
  testthat::expect_equal("list", class(draws))
  testthat::expect_equal("data.table", unique(sapply(draws, class)[1, ]))

  summarize_cols <- c("chain", "chain_draw", "draw")
  if (software == "tmb") summarize_cols <- c("draw")
  summary <- testthat::expect_silent(
    popMethods::popReconstruct_summarize_draws(
      draws = draws,
      summarize_cols = summarize_cols
    )
  )
  testthat::expect_equal("list", class(summary))
  testthat::expect_equal("data.table", unique(sapply(summary, class)[1, ]))
}

test_prior <- function(inputs, hyperparameters, settings, count_parameters) {
  draws <- testthat::expect_output(
    popMethods::popReconstruct_prior_draws(
      inputs = inputs,
      hyperparameters = hyperparameters,
      settings = settings,
      value_col = "value",
      method_name = "Original"
    )
  )
  testthat::expect_equal("list", class(draws))
  testthat::expect_equal("data.table", unique(sapply(draws, class)[1, ]))

  testthat::expect_silent(
    draws <- popMethods::popReconstruct_count_space_parameters(
      draws = draws,
      settings = settings,
      parameters = count_parameters,
      value_col = "value",
      quiet = TRUE
    )
  )
  testthat::expect_true(all(count_parameters %in% names(draws)))
  testthat::expect_equal("list", class(draws))
  testthat::expect_equal("data.table", unique(sapply(draws, class)[1, ]))

  summary_prior <- testthat::expect_silent(
    popMethods::popReconstruct_summarize_draws(
      draws = draws,
      summarize_cols = "draw"
    )
  )
  testthat::expect_equal("list", class(summary_prior))
  testthat::expect_equal("data.table", unique(sapply(summary_prior, class)[1, ]))
}

test_ccmpp_function <- function(inputs, settings) {

  population_demCore <- demCore::ccmpp(
    inputs = inputs,
    settings = settings
  )
  population_demCore <- demCore:::dt_to_matrix(
    dt = population_demCore,
    id_cols = c("year", "sex", "age_start", "age_end")
  )
  names(population_demCore) <- NULL
  for (i in 1:length(population_demCore)) {
    dimnames(population_demCore[[i]]) <- NULL

  }

  # convert from data.tables to matrices and list of matrices
  inputs_matrix <- lapply(inputs, function(dt) {
    input_matrix <- demCore:::dt_to_matrix(
      dt = dt,
      id_cols = setdiff(names(dt), "value")
    )
    return(input_matrix)
  })
  population_stan <- ccmpp(
    srb = inputs_matrix$srb,
    asfr = inputs_matrix$asfr,
    baseline = inputs_matrix$baseline,
    survival = inputs_matrix$survival,
    net_migration = inputs_matrix$net_migration,
    sexes = length(settings$sexes),
    interval = unique(diff(settings$ages)),
    A = length(settings$ages),
    A_m = length(settings$ages_mortality),
    A_f = length(settings$ages_asfr),
    A_f_offset = sum(settings$ages < min(settings$ages_asfr)),
    Y = length(settings$years)
  )

  testthat::expect_equal(population_demCore, population_stan)
}

test_nSx_function <- function(lt, settings) {

  id_cols <- c("year_start", "year_end", "sex", "age_start", "age_end")

  survival_demCore <- demCore::nSx_from_lx_nLx_Tx(
    dt = lt,
    id_cols = id_cols,
    terminal_age = max(settings$ages)
  )
  survival_demCore <- demCore:::dt_to_matrix(
    dt = survival_demCore[, .SD, .SDcols = c(id_cols, "nSx")],
    id_cols = id_cols, value_col = "nSx"
  )$female
  dimnames(survival_demCore) <- NULL

  survival_stan <- calculate_nSx(
    mx = demCore:::dt_to_matrix(
      dt = lt[, .SD, .SDcols = c(id_cols, "mx")],
      id_cols = id_cols, value_col = "mx"
    )$female,
    ax = demCore:::dt_to_matrix(
      dt = lt[, .SD, .SDcols = c(id_cols, "ax")],
      id_cols = id_cols, value_col = "ax"
    )$female,
    interval = unique(diff(settings$ages)),
    A_m = length(settings$ages_mortality),
    Y = length(settings$years)
  )

  testthat::expect_equal(survival_demCore, survival_stan)
}

# Test the original popReconstruct model specification --------------------

test_combinations(
  desc = "the original popReconstruct model works",
  inputs = demCore::burkina_faso_initial_estimates,
  data = demCore::burkina_faso_data,
  hyperparameters = hyperparameters,
  settings = settings,
  count_parameters = c("live_births", "net_migrants")
)

# Test stan function for ccmpp --------------------------------------------

description <- "ccmpp in stan gives the same output as the equivalent demCore function"
testthat::test_that("ccmpp in stan gives the same output as the equivalent demCore function", {
  test_ccmpp_function(
    inputs = demCore::burkina_faso_initial_estimates,
    settings = settings
  )
})

# run same tests but with ages_mortality = ages
new_settings <- copy(settings)
new_settings[["ages_mortality"]] <- new_settings[["ages"]]

new_inputs <- copy(demCore::burkina_faso_initial_estimates)
new_inputs[["survival"]] <- new_inputs[["survival"]][age_start <= max(new_settings$ages)]
new_inputs[["survival"]] <- new_inputs[["survival"]][
  age_start == max(new_settings$ages), age_end := Inf
]

description <- "ccmpp in stan gives the same output as the equivalent demCore function when 'ages_mortality=ages'"
testthat::test_that(description, {
  test_ccmpp_function(
    inputs = new_inputs,
    settings = new_settings
  )
})

# Test popReconstruct with aggregate population data ----------------------

# aggregate to total population
age_mapping <- data.table(age_start = 0, age_end = Inf)
aggregated_data <- copy(demCore::burkina_faso_data)
aggregated_data$population <- hierarchyUtils::agg(
  dt = aggregated_data$population,
  id_cols = c("year", "sex", "age_start", "age_end"),
  value_cols = "value",
  col_stem = "age",
  col_type = "interval",
  mapping = age_mapping
)

test_combinations(
  desc = "the popReconstruct model with aggregate data works",
  inputs = demCore::burkina_faso_initial_estimates,
  data = aggregated_data,
  hyperparameters = hyperparameters,
  settings = settings,
  count_parameters = c("live_births", "net_migrants"),
  test_software = c("stan", "tmb")
)

# Test popReconstruct with immigration/emigration -------------------------

# create inputs for immigration and emigration
new_inputs <- copy(demCore::burkina_faso_initial_estimates)
new_inputs$immigration <- new_inputs$net_migration
new_inputs$immigration[, value := 0.1]
new_inputs$emigration <- new_inputs$net_migration
new_inputs$emigration[, value := 0.1]
new_inputs$net_migration <- NULL

new_hyperparameters <- copy(hyperparameters)
new_hyperparameters$immigration <- hyperparameters$net_migration
new_hyperparameters$emigration <- hyperparameters$net_migration
new_hyperparameters$net_migration <- NULL

test_combinations(
  desc = "the popReconstruct model with immigration/emigration parameters works",
  inputs = new_inputs,
  data = demCore::burkina_faso_data,
  hyperparameters = new_hyperparameters,
  settings = settings,
  count_parameters = c("live_births", "immigrants", "emigrants")
)

# Test popReconstruct with mx, ax -----------------------------------------

# create rough inputs for mx and ax
lt <- copy(demCore::burkina_faso_initial_estimates$survival)
# px approximately equal to survivorship ratio
lt[, qx := 1 - value]
lt[, ax := 2.5]
lt[, value := NULL]
lt[is.infinite(age_end), c("qx", "ax") := list(1, 5)]
id_cols <- c("year_start", "year_end", "sex", "age_start", "age_end")
demCore::lifetable(lt, id_cols = id_cols)
lt[, age_length := NULL]

# add inputs for mx and ax
new_inputs <- copy(demCore::burkina_faso_initial_estimates)
new_inputs$mx <- lt[, .SD, .SDcols = c(id_cols, "mx")]
setnames(new_inputs$mx, "mx", "value")
new_inputs$ax <- lt[, .SD, .SDcols = c(id_cols, "ax")]
setnames(new_inputs$ax, "ax", "value")
new_inputs$survival <- NULL

hyperparameters_mx_ax <- copy(hyperparameters)
hyperparameters_mx_ax$mx <- hyperparameters$survival
hyperparameters_mx_ax$survival <- NULL
hyperparameters_mx_ax$mx$beta <- 0.000109
hyperparameters_mx_ax$non_terminal_ax <- hyperparameters_mx_ax$mx
hyperparameters_mx_ax$terminal_ax <- hyperparameters_mx_ax$mx

test_combinations(
  desc = "the popReconstruct model with mx/ax parameters works",
  inputs = new_inputs,
  data = demCore::burkina_faso_data,
  hyperparameters = hyperparameters_mx_ax,
  settings = settings,
  count_parameters = c("live_births", "deaths", "net_migrants")
)

test_nSx_function(lt, settings)

# Test popReconstruct with mx ---------------------------------------------

hyperparameters_mx <- copy(hyperparameters)
hyperparameters_mx$mx <- hyperparameters$survival
hyperparameters_mx$mx$beta <- 0.000109
hyperparameters_mx$survival <- NULL

# don't estimate ax
settings_mx <- copy(settings)
settings_mx$fixed_parameters <- c("non_terminal_ax", "terminal_ax")

test_combinations(
  desc = "the popReconstruct model with mx parameter works",
  inputs = new_inputs,
  data = demCore::burkina_faso_data,
  hyperparameters = hyperparameters_mx,
  settings = settings_mx,
  count_parameters = c("live_births", "deaths", "net_migrants")
)

# Test popReconstruct with mx, ax (ages_mortality = ages) -----------------

# collapse life table inputs to have the same age groups as other parameters
mapping <- data.table(age_start = settings$ages)
hierarchyUtils::gen_end(mapping, id_cols = "age_start", col_stem = "age")
lt <- demCore::agg_lt(
  dt = lt, id_cols = id_cols,
  age_mapping = mapping, quiet = TRUE,
  present_agg_severity = "none"
)
hierarchyUtils::gen_length(lt, col_stem = "age")
lt[, mx := demCore::qx_ax_to_mx(qx, ax, age_length)]

# add inputs for mx and ax
new_inputs <- copy(demCore::burkina_faso_initial_estimates)
new_inputs$mx <- lt[, .SD, .SDcols = c(id_cols, "mx")]
setnames(new_inputs$mx, "mx", "value")
new_inputs$ax <- lt[, .SD, .SDcols = c(id_cols, "ax")]
setnames(new_inputs$ax, "ax", "value")
new_inputs$survival <- NULL

new_settings <- copy(settings)
new_settings$ages_mortality <- new_settings$ages

test_combinations(
  desc = "the popReconstruct model with mx/ax parameters & 'ages_mortality=ages' works",
  inputs = new_inputs,
  data = demCore::burkina_faso_data,
  hyperparameters = hyperparameters_mx_ax,
  settings = new_settings,
  count_parameters = c("live_births", "deaths", "net_migrants")
)

demCore::lifetable(lt, id_cols)
test_nSx_function(lt, new_settings)
