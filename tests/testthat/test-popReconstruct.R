settings = list(
  years = seq(1960, 2000, 5),
  sexes = c("female"),
  ages = seq(0, 80, 5),
  ages_survival = seq(0, 85, 5),
  ages_asfr = seq(15, 45, 5),

  n_draws = 100
)
hyperparameters <- list(
  asfr = list(alpha = 1, beta = 0.0109),
  population = list(alpha = 1, beta = 0.0109),
  survival = list(alpha = 1, beta = 0.0109),
  net_migration = list(alpha = 1, beta = 0.0436)
)


test_fit <- function(inputs, hyperparameters, software, ...) {

  fit <- testthat::expect_error(
    popMethods::popReconstruct_fit(
      inputs = inputs,
      data = demCore::burkina_faso_data,
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

test_prior <- function(inputs, hyperparameters) {
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
    hyperparameters = hyperparameters,
    software = "stan",
    chains = 1, warmup = 100, iter = 200, thin = 2, seed = 3
  )
})

testthat::test_that("the original popReconstruct model works in tmb", {
  test_fit(
    inputs = demCore::burkina_faso_initial_estimates,
    hyperparameters = hyperparameters,
    software = "tmb"
  )
})

testthat::test_that("sampling from the original popReconstruct model prior works", {
  test_prior(
    inputs = demCore::burkina_faso_initial_estimates,
    hyperparameters = hyperparameters
  )
})

# Test popReconstruct with immigration/emigration -------------------------

# create inputs for immigration and emigration
inputs <- copy(demCore::burkina_faso_initial_estimates)
inputs$immigration <- inputs$net_migration
inputs$immigration[, value := 0.1]
inputs$emigration <- inputs$net_migration
inputs$emigration[, value := 0.1]
inputs$net_migration <- NULL

hyperparameters$immigration <- hyperparameters$net_migration
hyperparameters$emigration <- hyperparameters$net_migration
hyperparameters$net_migration <- NULL

testthat::test_that("the popReconstruct (immigration/emigration)  model works in stan", {
  test_fit(
    inputs = inputs,
    hyperparameters = hyperparameters,
    software = "stan",
    chains = 1, warmup = 100, iter = 200, thin = 2, seed = 3
  )
})

testthat::test_that("the popReconstruct model (immigration/emigration) works in tmb", {
  test_fit(
    inputs = inputs,
    hyperparameters = hyperparameters,
    software = "tmb"
  )
})

testthat::test_that("sampling from popReconstruct (immigration/emigration) model prior works", {
  test_prior(
    inputs = inputs,
    hyperparameters = hyperparameters
  )
})
