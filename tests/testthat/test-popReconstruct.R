

# Test the original popReconstruct model specification --------------------

original_settings = list(
  years = seq(1960, 2000, 5),
  sexes = c("female"),
  ages = seq(0, 80, 5),
  ages_survival = seq(0, 85, 5),
  ages_asfr = seq(15, 45, 5),

  n_draws = 100
)
original_hyperparameters <- list(asfr = list(alpha = 1, beta = 0.0109),
                                 population = list(alpha = 1, beta = 0.0109),
                                 survival = list(alpha = 1, beta = 0.0109),
                                 net_migration = list(alpha = 1, beta = 0.0436))

testthat::test_that("the original popReconstruct model works in stan", {
  fit_stan <- testthat::expect_error(
    popMethods::popReconstruct_fit(
      inputs = demCore::burkina_faso_initial_estimates,
      data = demCore::burkina_faso_data,
      hyperparameters = original_hyperparameters,
      settings = original_settings,
      value_col = "value",
      software = "stan",
      chains = 2, warmup = 100, iter = 200, thin = 2, seed = 3
    ),
    NA
  )

  draws_stan <- testthat::expect_silent(
    popMethods::popReconstruct_posterior_draws(
      fit = fit_stan,
      inputs = demCore::burkina_faso_initial_estimates,
      settings = original_settings,
      value_col = "value",
      software = "stan",
      method_name = "original"
    )
  )
  testthat::expect_equal("list", class(draws_stan))
  testthat::expect_equal("data.table", unique(sapply(draws_stan, class)[1, ]))

  summary_stan <- testthat::expect_silent(
    popMethods::popReconstruct_summarize_draws(
      draws = draws_stan
    )
  )
  testthat::expect_equal("list", class(summary_stan))
  testthat::expect_equal("data.table", unique(sapply(summary_stan, class)[1, ]))

})

testthat::test_that("the original popReconstruct model works in tmb", {
  fit_tmb <- testthat::expect_error(
    popMethods::popReconstruct_fit(
      inputs = demCore::burkina_faso_initial_estimates,
      data = demCore::burkina_faso_data,
      hyperparameters = original_hyperparameters,
      settings = original_settings,
      value_col = "value",
      software = "tmb",
    ),
    NA
  )

  draws_tmb <- testthat::expect_silent(
    popMethods::popReconstruct_posterior_draws(
      fit = fit_tmb,
      inputs = demCore::burkina_faso_initial_estimates,
      settings = original_settings,
      value_col = "value",
      software = "tmb",
      method_name = "original"
    )
  )
  testthat::expect_equal("list", class(draws_tmb))
  testthat::expect_equal("data.table", unique(sapply(draws_tmb, class)[1, ]))

  summary_tmb <- testthat::expect_silent(
    popMethods::popReconstruct_summarize_draws(
      draws = draws_tmb,
      summarize_cols = "draw"
    )
  )
  testthat::expect_equal("list", class(summary_tmb))
  testthat::expect_equal("data.table", unique(sapply(summary_tmb, class)[1, ]))
})

testthat::test_that("sampling from the original popReconstruct model prior works", {
  draws_prior <- testthat::expect_output(
    popMethods::popReconstruct_prior_draws(
      inputs = demCore::burkina_faso_initial_estimates,
      hyperparameters = original_hyperparameters,
      settings = original_settings,
      value_col = "value",
      method_name = "Original"
    )
  )
  testthat::expect_equal("list", class(draws_prior))
  testthat::expect_equal("data.table", unique(sapply(draws_prior, class)[1, ]))

  summary_prior <- testthat::expect_silent(
    popMethods::popReconstruct_summarize_draws(
      draws = draws_prior,
      summarize_cols = "draw"
    )
  )
  testthat::expect_equal("list", class(summary_prior))
  testthat::expect_equal("data.table", unique(sapply(summary_prior, class)[1, ]))
})
