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

test_fit <- function(inputs, data, hyperparameters, settings, software, ...) {

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

test_prior <- function(inputs, hyperparameters, settings) {
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

  summary_prior <- testthat::expect_silent(
    popMethods::popReconstruct_summarize_draws(
      draws = draws,
      summarize_cols = "draw"
    )
  )
  testthat::expect_equal("list", class(summary_prior))
  testthat::expect_equal("data.table", unique(sapply(summary_prior, class)[1, ]))
}

# Test the original popReconstruct model specification --------------------

testthat::test_that("the original popReconstruct model works in stan", {
  test_fit(
    inputs = demCore::burkina_faso_initial_estimates,
    data = demCore::burkina_faso_data,
    hyperparameters = hyperparameters,
    settings = settings,
    software = "stan",
    chains = 1, warmup = 100, iter = 200, thin = 2, seed = 3
  )
})

testthat::test_that("the original popReconstruct model works in tmb", {
  test_fit(
    inputs = demCore::burkina_faso_initial_estimates,
    data = demCore::burkina_faso_data,
    hyperparameters = hyperparameters,
    settings = settings,
    software = "tmb"
  )
})

testthat::test_that("sampling from the original popReconstruct model prior works", {
  test_prior(
    inputs = demCore::burkina_faso_initial_estimates,
    hyperparameters = hyperparameters,
    settings = settings
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

testthat::test_that("the popReconstruct model (with aggregate data) works in stan", {
  test_fit(
    inputs = demCore::burkina_faso_initial_estimates,
    data = aggregated_data,
    hyperparameters = hyperparameters,
    settings = settings,
    software = "stan",
    chains = 1, warmup = 100, iter = 200, thin = 2, seed = 3
  )
})

testthat::test_that("the popReconstruct model (with aggregate data) works in tmb", {
  test_fit(
    inputs = demCore::burkina_faso_initial_estimates,
    data = aggregated_data,
    hyperparameters = hyperparameters,
    settings = settings,
    software = "tmb"
  )
})

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

testthat::test_that("the popReconstruct (immigration/emigration)  model works in stan", {
  test_fit(
    inputs = new_inputs,
    data = demCore::burkina_faso_data,
    hyperparameters = new_hyperparameters,
    settings = settings,
    software = "stan",
    chains = 1, warmup = 100, iter = 200, thin = 2, seed = 3
  )
})

testthat::test_that("the popReconstruct model (immigration/emigration) works in tmb", {
  test_fit(
    inputs = new_inputs,
    data = demCore::burkina_faso_data,
    hyperparameters = new_hyperparameters,
    settings = settings,
    software = "tmb"
  )
})

testthat::test_that("sampling from popReconstruct (immigration/emigration) model prior works", {
  test_prior(
    inputs = new_inputs,
    hyperparameters = new_hyperparameters,
    settings = settings
  )
})

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


testthat::test_that("the popReconstruct (mx & ax) model works in stan", {
  test_fit(
    inputs = new_inputs,
    data = demCore::burkina_faso_data,
    hyperparameters = hyperparameters_mx_ax,
    settings = settings,
    software = "stan",
    chains = 1, warmup = 100, iter = 200, thin = 2, seed = 3
  )
})

testthat::test_that("the popReconstruct model (mx & ax) works in tmb", {
  test_fit(
    inputs = new_inputs,
    data = demCore::burkina_faso_data,
    hyperparameters = hyperparameters_mx_ax,
    settings = settings,
    software = "tmb"
  )
})

testthat::test_that("sampling from popReconstruct (mx & ax) model prior works", {
  test_prior(
    inputs = new_inputs,
    hyperparameters = hyperparameters_mx_ax,
    settings = settings
  )
})

# Test popReconstruct with mx ---------------------------------------------

hyperparameters_mx <- copy(hyperparameters)
hyperparameters_mx$mx <- hyperparameters$survival
hyperparameters_mx$mx$beta <- 0.000109
hyperparameters_mx$survival <- NULL

# don't estimate ax
settings_mx <- copy(settings)
settings_mx$fixed_parameters <- c("non_terminal_ax", "terminal_ax")

testthat::test_that("the popReconstruct (mx) model works in stan", {
  test_fit(
    inputs = new_inputs,
    data = demCore::burkina_faso_data,
    hyperparameters = hyperparameters_mx,
    settings = settings_mx,
    software = "stan",
    chains = 1, warmup = 100, iter = 200, thin = 2, seed = 3
  )
})

testthat::test_that("the popReconstruct model (mx) works in tmb", {
  test_fit(
    inputs = new_inputs,
    data = demCore::burkina_faso_data,
    hyperparameters = hyperparameters_mx,
    settings = settings_mx,
    software = "tmb"
  )
})

testthat::test_that("sampling from popReconstruct (mx) model prior works", {
  test_prior(
    inputs = new_inputs,
    hyperparameters = hyperparameters_mx,
    settings = settings_mx
  )
})
