#' @export
#' @rdname popReconstruct_fit
#' @include popReconstruct_fit.R
popReconstruct_posterior_draws <- function(fit,
                                           inputs,
                                           settings,
                                           value_col,
                                           software,
                                           method_name) {

  # Validate arguments ------------------------------------------------------

  # check `software` argument
  assertthat::assert_that(
    assertthat::is.string(software),
    software %in% c("stan", "tmb"),
    msg = "`software` must be 'stan' or 'tmb'."
  )

  # check `method_name` argument
  assertthat::assert_that(
    assertthat::is.string(method_name),
    msg = "`method_name` must be a character vector of length one."
  )

  demCore:::validate_ccmpp_inputs(inputs, settings, value_col)

  # add optional settings needed for popReconstruct
  settings <- copy(settings)
  settings <- create_optional_settings(settings, inputs)
  detailed_settings <- create_detailed_settings(settings)

  validate_popReconstruct_inputs(inputs, settings, detailed_settings, value_col)

  # Extract modeled parameters ----------------------------------------------

  # extract directly modeled offset parameters
  if (software == "stan") {
    assertthat::assert_that(
      class(fit) == "stanfit",
      msg = "'fit' object is not a stan object"
    )
    draws <- extract_stan_draws(fit, inputs, settings, detailed_settings)
  } else if (software == "tmb") {
    assertthat::assert_that(
      class(fit) == "sdreport",
      msg = "'fit' object is not a tmb object"
    )
    draws <- extract_tmb_draws(fit, inputs, value_col, settings, detailed_settings)
  }

  # add method column
  for (measure in names(draws)) {
    draws[[measure]][, method := method_name]
    data.table::setnames(draws[[measure]], "value", value_col)
  }

  return(draws)
}

#' @title Helper function to extract draws from popReconstruct stan fits
#'
#' @inheritParams popReconstruct_posterior_draws
#' @param detailed_settings \[`list()`\]\cr
#'   Detailed settings for each ccmpp input and 'population'.
#'
#' @seealso [`popReconstruct_posterior_draws()`]
extract_stan_draws <- function(fit, inputs, settings, detailed_settings) {

  # compile draws from each of the chains
  draws <- rstan::extract(fit, permuted = FALSE)
  chains <- dim(draws)[2]
  draws <- lapply(1:chains, function(chain) {
    temp <- data.table(draws[, chain, ])
    temp[, chain := as.integer(chain)]
    temp[, chain_draw := .I]
    temp <- melt(temp, id.vars = c("chain", "chain_draw"),
                 variable.name = "parameter",
                 value.name = "value")
    return(temp)
  })
  draws <- rbindlist(draws)
  draws[, draw := as.integer(((chain - 1) * max(chain_draw)) + chain_draw)]

  format_draws <- function(comp, param, settings, comp_detailed_settings) {

    # get component specific draws
    measure_type <- comp_detailed_settings[["measure_type"]]
    id_cols <- comp_detailed_settings[["id_cols"]]
    if (grepl("^offset", param) & comp != "baseline") {
      years <- comp_detailed_settings[["years_knots"]]
    } else {
      years <- comp_detailed_settings[["years"]]
    }
    years_projections <- settings[["years_projections"]]
    int <- settings[["int"]]
    sexes <- comp_detailed_settings[["sexes"]]
    if (grepl("^offset", param)) {
      ages <- comp_detailed_settings[["ages_knots"]]
    } else {
      ages <- comp_detailed_settings[["ages"]]
    }

    # subset to just the draws for this specific parameter
    component_draws <- draws[grepl(paste0("^", param), parameter)]

    # format the draws if they were actually estimated
    if (nrow(component_draws) > 0) {
      # add on id variables to draws
      component_draws[, parameter := gsub(paste0("^", param, "|\\[|\\]"), "", parameter)]
      component_draws[, c(if (grepl("^offset", param) | param %in% c("immigration", "emigration")) "estimate",
                          if (!is.null(sexes)) "sex_index", "age_index", "year_index") :=
                        data.table::tstrsplit(parameter, split = ",")]
      assertthat::assert_that(
        assertthat::are_equal(
          length(unique(component_draws$year_index)),
          length(years)
        )
      )
      component_draws[, year_start := years[as.integer(year_index)]]
      if (!is.null(sexes)) {
        assertthat::assert_that(
          assertthat::are_equal(
            length(unique(component_draws$sex_index)),
            length(sexes)
          )
        )
        component_draws[, sex := sexes[as.integer(sex_index)]]
      }
      if (!is.null(ages)) {
        assertthat::assert_that(
          assertthat::are_equal(
            length(unique(component_draws$age_index)),
            length(ages)
          )
        )
        component_draws[, age_start := ages[as.integer(age_index)]]
      }

    } else {
      # create draws of zeroes for offset parameters that were fixed
      component_draws <- list(years, sexes, ages,
                              chain = unique(draws$chain),
                              chain_draw = unique(draws$chain_draw),
                              value = 0)
      names(component_draws)[1:3] <- c("year_start", "sex", "age_start")
      component_draws <- component_draws[!mapply(is.null, component_draws)]
      component_draws <- do.call(data.table::CJ, component_draws)
      component_draws[, draw := ((chain - 1) * max(chain_draw)) + chain_draw]
    }

    if (measure_type == "stock") setnames(component_draws, "year_start", "year")

    non_end_id_cols <- c(id_cols[!grepl("_end$", id_cols)],
                         "chain", "chain_draw", "draw")
    component_draws <- component_draws[, c(non_end_id_cols, "value"), with = F]

    # add on the 'year_end' column
    if ("year_start" %in% id_cols) {
      hierarchyUtils::gen_end(
        dt = component_draws,
        id_cols = non_end_id_cols,
        col_stem = "year",
        right_most_endpoint = max(years_projections)
      )
    }
    # add on the 'age_end' column
    if ("age_start" %in% id_cols) {
      hierarchyUtils::gen_end(
        dt = component_draws,
        id_cols = non_end_id_cols,
        col_stem = "age",
        right_most_endpoint = ifelse(comp == "asfr", max(ages) + int, Inf)
      )
    }

    setkeyv(component_draws, id_cols)
    return(component_draws)
  }

  # extract offset parameters
  offset_draws <- lapply(names(inputs), function(comp) {
    transformation_name <- detailed_settings[[comp]][["transformation_name"]]
    comp_offset_draws <- format_draws(
      comp = comp,
      param = paste0("offset_",
                     if (!is.null(transformation_name)) transformation_name,
                     if (!is.null(transformation_name)) "_",
                     comp),
      settings,
      detailed_settings[[comp]]
    )
    return(comp_offset_draws)
  })
  names(offset_draws) <- paste0("offset_", names(inputs))

  # extract spline offset parameters
  spline_offset_draws <- lapply(names(inputs), function(comp) {
    transformation_name <- detailed_settings[[comp]][["transformation_name"]]
    comp_spline_offset_draws <- format_draws(
      comp = comp,
      param = paste0("spline_offset_",
                     if (!is.null(transformation_name)) transformation_name,
                     if (!is.null(transformation_name)) "_",
                     comp),
      settings,
      detailed_settings[[comp]]
    )
    return(comp_spline_offset_draws)
  })
  names(spline_offset_draws) <- paste0("spline_offset_", names(inputs))

  # extract ccmpp input parameters
  ccmpp_input_draws <- lapply(names(inputs), function(comp) {
    comp_draws <- format_draws(
      comp = comp,
      param = paste0(comp),
      settings,
      detailed_settings[[comp]]
    )
    return(comp_draws)
  })
  names(ccmpp_input_draws) <- names(inputs)

  # extract projected population parameters
  population_draws <- format_draws(
    comp = "population",
    param = "population\\[",
    settings,
    detailed_settings[["population"]]
  )

  # extract variance parameters
  variance_draws <- draws[grepl("^sigma2", parameter)]
  variance_draws[, parameter := gsub("^sigma2_|\\[1\\]", "", parameter)]
  variance_draws <- variance_draws[, list(parameter, chain, chain_draw,
                                          draw, value)]
  setkeyv(variance_draws, c("parameter", "chain", "chain_draw", "draw"))

  draws <- list(variance = variance_draws, population = population_draws)
  draws <- c(draws, offset_draws, spline_offset_draws, ccmpp_input_draws)
  return(draws)
}

#' @title Extract draws from the popReconstruct model TMB fit
#'
#' @inheritParams popReconstruct_posterior_draws
#' @inheritParams extract_stan_draws
#'
#' @seealso [`popReconstruct_posterior_draws()`]
extract_tmb_draws <- function(fit, inputs, value_col, settings, detailed_settings) {

  # Generate draws for random and fixed parameters --------------------------

  # draw from a multivariate normal distribution given the mean and precision
  # matrix for all parameters
  gen_draws <- function(mu, prec, n_draws) {
    z = matrix(stats::rnorm(length(mu) * n_draws), ncol = n_draws)
    L_inv = Matrix::Cholesky(prec)
    mu + solve(as(L_inv, "pMatrix"), solve(Matrix::t(as(L_inv, "Matrix")), z))
  }

  # extract mean and precision matrix for all model parameters
  random_mean <- fit$par.random
  ii <- grep("^offset_", rownames(fit$jointPrecision))
  random_prec <- fit$jointPrecision[ii,ii]
  random_draws <- gen_draws(
    mu = random_mean,
    prec = random_prec,
    n_draws = settings$n_draws
  )

  fixed_mean <- fit$par.fixed
  ii <- grep("^log_sigma", rownames(fit$jointPrecision))
  fixed_prec <- fit$jointPrecision[ii,ii]
  fixed_draws <- gen_draws(
    mu = fixed_mean,
    prec = fixed_prec,
    n_draws = settings$n_draws
  )

  # Format the offset and variance draws ------------------------------------

  # format variance draws
  variance_draws <- data.table(as.matrix(fixed_draws))
  variance_draws[, parameter := gsub("log_sigma2_", "", names(fixed_mean))]
  variance_draws <- melt(
    variance_draws,
    id.vars = "parameter",
    variable.name = "draw",
    value.name = "value"
  )
  variance_draws[, draw := as.integer(draw)]
  variance_draws[, value := exp(value)]
  setkeyv(variance_draws, c("parameter", "draw"))

  format_draws <- function(comp, param, settings, comp_detailed_settings) {

    # get component specific draws
    measure_type <- comp_detailed_settings[["measure_type"]]
    id_cols <- comp_detailed_settings[["id_cols"]]
    if (grepl("^offset", param) & comp != "baseline") {
      years <- comp_detailed_settings[["years_knots"]]
    } else {
      years <- comp_detailed_settings[["years"]]
    }
    years_projections <- settings[["years_projections"]]
    int <- settings[["int"]]
    sexes <- comp_detailed_settings[["sexes"]]
    if (grepl("^offset", param)) {
      ages <- comp_detailed_settings[["ages_knots"]]
    } else {
      ages <- comp_detailed_settings[["ages"]]
    }

    # subset to just the draws for this specific parameter
    component_draws <- as.data.table(random_draws[param == names(random_mean),])

    # format the draws if they were actually estimated
    if (nrow(component_draws) > 0) {
      # add on id variables to draws
      component_id_vars <- list(year_start = years, sex = sexes,
                                age_start = ages)
      component_id_vars <- component_id_vars[!mapply(is.null, component_id_vars)]
      component_id_vars <- do.call(data.table::CJ, component_id_vars)
      component_draws <- cbind(component_id_vars, component_draws)

      # switch to long format
      component_draws <- melt(
        component_draws,
        id.vars = names(component_id_vars),
        variable.name = "draw",
        value.name = "value"
      )
      component_draws[, draw := as.integer(draw)]

    } else {
      # create draws of zeroes for offset parameters that were fixed
      component_draws <- list(year_start = years, sex = sexes, age_start = ages,
                              draw = 1:settings$n_draws,
                              value = 0)
      names(component_draws)[1:3] <- c("year_start", "sex", "age_start")
      component_draws <- component_draws[!mapply(is.null, component_draws)]
      component_draws <- do.call(data.table::CJ, component_draws)
    }

    if (measure_type == "stock") setnames(component_draws, "year_start", "year")

    non_end_id_cols <- c(id_cols[!grepl("_end$", id_cols)], "draw")
    component_draws <- component_draws[, c(non_end_id_cols, "value"), with = F]

    # add on the 'year_end' column
    if ("year_start" %in% id_cols) {
      hierarchyUtils::gen_end(
        dt = component_draws,
        id_cols = non_end_id_cols,
        col_stem = "year",
        right_most_endpoint = max(settings$years_projections)
      )
    }
    # add on the 'age_end' column
    if ("age_start" %in% id_cols) {
      hierarchyUtils::gen_end(
        dt = component_draws,
        id_cols = non_end_id_cols,
        col_stem = "age",
        right_most_endpoint = ifelse(comp == "asfr", max(ages) + int, Inf)
      )
    }

    setkeyv(component_draws, id_cols)
    return(component_draws)
  }

  # extract offset parameters
  offset_draws <- lapply(names(inputs), function(comp) {
    transformation_name <- detailed_settings[[comp]][["transformation_name"]]
    comp_offset_draws <- format_draws(
      comp = comp,
      param = paste0("offset_",
                     if (!is.null(transformation_name)) transformation_name,
                     if (!is.null(transformation_name)) "_",
                     comp),
      settings,
      detailed_settings[[comp]]
    )
    return(comp_offset_draws)
  })
  names(offset_draws) <- names(inputs)

  # Calculate derived parameters --------------------------------------------

  spline_offset_draws <- calculate_spline_offset_draws(
    offset_draws,
    detailed_settings,
    draw_col_name = "draw"
  )

  ccmpp_input_draws <- calculate_ccmpp_input_draws(
    spline_offset_draws,
    inputs,
    detailed_settings,
    value_col
  )

  # project population using ccmpp at draw level
  pop_draws <- ccmpp_draws(
    ccmpp_input_draws,
    settings,
    settings$n_draws,
    draw_col_name = "draw"
  )
  # create 'age_end' column
  hierarchyUtils::gen_end(
    dt = pop_draws,
    id_cols = c("year", "sex", "age_start", "draw"),
    col_stem = "age"
  )

  # Combine results ---------------------------------------------------------

  draws <- list(variance = variance_draws,
                population = pop_draws)
  names(offset_draws) <- paste0("offset_", names(offset_draws))
  names(spline_offset_draws) <- paste0("spline_offset_",
                                       names(spline_offset_draws))
  draws <- c(draws, ccmpp_input_draws, offset_draws, spline_offset_draws)
  return(draws)
}

#' @export
#' @rdname popReconstruct_fit
popReconstruct_prior_draws <- function(inputs,
                                       hyperparameters,
                                       settings,
                                       value_col,
                                       method_name,
                                       chunk_size = 100) {

  # Validate arguments ------------------------------------------------------

  demCore:::validate_ccmpp_inputs(inputs, settings, value_col)

  # check `chunk_size` argument
  assertthat::assert_that(
    assertthat::is.count(chunk_size),
    chunk_size <= settings$n_draws,
    msg = "`chunk_size` must be an integer less than the 'n_draws' setting."
  )

  # check `method_name` argument
  assertthat::assert_that(
    assertthat::is.string(method_name),
    msg = "`method_name` must be a character vector of length one."
  )

  # add optional settings needed for popReconstruct
  settings <- copy(settings)
  settings <- create_optional_settings(settings, inputs)
  detailed_settings <- create_detailed_settings(settings)

  validate_popReconstruct_hyperparameters(hyperparameters, inputs, settings)
  validate_popReconstruct_inputs(inputs, settings, detailed_settings, value_col)

  # Sample from prior distribution ------------------------------------------

  print(paste("Sampling", settings$n_draws, "draws from the prior distribution",
              chunk_size, "draws at a time"))

  # create overall containers for draws
  all_variance <- NULL
  all_offset_draws <- NULL
  all_spline_offset_draws <- NULL
  all_ccmpp_input_draws <- NULL
  all_pop_draws <- NULL

  positive_draws_created <- 0
  while(positive_draws_created < settings$n_draws) {

    # sample parameter specific variance draws from inverse gamma distribution
    print("- Sampling variance draws")
    variance <- lapply(settings$estimated_parameters, function(parameter) {
      if (parameter == "baseline") parameter <- "population"
      draws <- data.table(parameter = parameter,
                          original_draw = 1:chunk_size,
                          value = demUtils::rinvgamma(
                            n = chunk_size,
                            shape = hyperparameters[[parameter]][["alpha"]],
                            scale = hyperparameters[[parameter]][["beta"]])
      )
    })
    variance <-data.table::rbindlist(variance)

    # sample ccmpp input offset parameters
    print("- Sampling offset parameter draws")
    offset_draws <- lapply(names(inputs), function(comp) {

      id_cols <- detailed_settings[[comp]][["id_cols"]]

      offset_component <- lapply(1:chunk_size, function(d) {

        # create id variables only data.table
        combinations <- list(
          if ("year_start" %in% id_cols) detailed_settings[[comp]][["years_knots"]],
          if (comp == "baseline") detailed_settings[[comp]][["years"]],
          if ("sex" %in% id_cols) settings$sexes,
          if ("age_start" %in% id_cols) detailed_settings[[comp]][["ages_knots"]],
          original_draw = d,
          value = 0
        )
        names(combinations)[1:4] <- c("year_start", "year", "sex", "age_start")
        combinations <- combinations[!mapply(is.null, combinations)]
        offset <- do.call(data.table::CJ, combinations)

        # sample the actual offset values if component is not fixed
        if (comp %in% settings$estimated_parameters) {
          var_comp <- ifelse(comp == "baseline", "population", comp)
          var_val <- variance[parameter == var_comp & original_draw == d, value]
          offset[, value := stats::rnorm(n = .N, mean = 0, sd = sqrt(var_val))]
        }
        return(offset)
      })
      offset_component <- data.table::rbindlist(offset_component)

      # add on the 'year_end' column
      if ("year_start" %in% id_cols) {
        hierarchyUtils::gen_end(
          dt = offset_component,
          id_cols = c(id_cols[!grepl("_end$", id_cols)], "original_draw"),
          col_stem = "year",
          right_most_endpoint = ifelse(comp == "baseline",
                                       min(settings$years),
                                       max(settings$years_projections))
        )
      }
      # add on the 'age_end' column
      if ("age_start" %in% id_cols) {
        hierarchyUtils::gen_end(
          dt = offset_component,
          id_cols = c(id_cols[!grepl("_end$", id_cols)], "original_draw"),
          col_stem = "age",
          right_most_endpoint = ifelse(comp == "asfr",
                                       max(settings$ages_asfr) + settings$int,
                                       Inf)
        )
      }
      return(offset_component)
    })
    names(offset_draws) <- names(inputs)

    # calculate spline offset draws
    print("- Calculating spline offset draws")
    spline_offset_draws <- calculate_spline_offset_draws(
      offset_draws,
      detailed_settings,
      draw_col_name = "original_draw"
    )

    print("- Calculating ccmpp input draws")
    ccmpp_input_draws <- calculate_ccmpp_input_draws(
      spline_offset_draws,
      inputs,
      detailed_settings,
      value_col
    )

    # project population using ccmpp at draw level
    print("- Projecting population draws using ccmpp")
    pop_draws <- ccmpp_draws(
      ccmpp_input_draws,
      settings,
      chunk_size,
      draw_col_name = "original_draw"
    )

    print("- Subsetting to draws that satisfy the positive population constraint")
    # determine which draws satisfy the constraint that the projected population for all years-sexes-ages is positive
    negative_draws <- pop_draws[value < 0, unique(original_draw)]
    positive_draws <- pop_draws[!original_draw %in% negative_draws, unique(original_draw)]
    if (length(positive_draws) == 0) next()

    # renumber draws starting from the last
    new_draw_numbers <- (positive_draws_created + 1):((positive_draws_created) + length(positive_draws))
    new_draw_numbers <- data.table(original_draw = positive_draws, draw = new_draw_numbers)

    # append draws that satisfy positive population constraint
    update_draw_numbers <- function(draws, positive_draws, new_draw_numbers) {
      if (any(class(draws) == "list")) {
        # determine the named elements corresponding to the estimate tables
        measures <- names(draws)

        # combine together all estimates from different methods
        combined_results <- lapply(measures, function(measure){
          measure_draws <- draws[[measure]][original_draw %in% positive_draws]
          measure_draws <- merge(measure_draws, new_draw_numbers, by = "original_draw", all = T)
          measure_draws[, original_draw := NULL]
          return(measure_draws)
        })
        names(combined_results) <- measures
      } else {
        combined_results <- draws[original_draw %in% positive_draws]
        combined_results <- merge(combined_results, new_draw_numbers, by = "original_draw", all = T)
        combined_results[, original_draw := NULL]
      }
      return(combined_results)
    }
    all_variance <- rbind(
      all_variance,
      update_draw_numbers(variance, positive_draws, new_draw_numbers)
    )
    all_offset_draws <- rbindlist_dts(
      all_offset_draws,
      update_draw_numbers(offset_draws, positive_draws, new_draw_numbers)
    )
    all_spline_offset_draws <- rbindlist_dts(
      all_spline_offset_draws,
      update_draw_numbers(spline_offset_draws, positive_draws, new_draw_numbers)
    )
    all_ccmpp_input_draws <- rbindlist_dts(
      all_ccmpp_input_draws,
      update_draw_numbers(ccmpp_input_draws, positive_draws, new_draw_numbers)
    )
    all_pop_draws <- rbind(
      all_pop_draws,
      update_draw_numbers(pop_draws, positive_draws, new_draw_numbers)
    )

    positive_draws_created <- positive_draws_created + length(positive_draws)
    print(paste(min(positive_draws_created, settings$n_draws), "out of",
                settings$n_draws, "total draws that satisfy the positive",
                "population constraint created"))
  }

  # create 'age_end' column
  hierarchyUtils::gen_end(
    dt = all_pop_draws,
    id_cols = c("year", "sex", "age_start", "draw"),
    col_stem = "age"
  )

  # combine all results together
  draws <- list(variance = all_variance,
                population = all_pop_draws)
  names(all_offset_draws) <- paste0("offset_", names(all_offset_draws))
  names(all_spline_offset_draws) <- paste0("spline_offset_",
                                           names(all_spline_offset_draws))
  draws <- c(draws, all_offset_draws, all_spline_offset_draws,
             all_ccmpp_input_draws)

  # drop extra draws that may have been sampled
  # add method column
  for (measure in names(draws)) {
    draws[[measure]] <- draws[[measure]][draw <= settings$n_draws]
    draws[[measure]][, method := method_name]
    data.table::setnames(draws[[measure]], "value", value_col)
  }

  return(draws)
}

#' @title Helper function to calculate spline offsets for all ccmpp inputs
#'
#' @param offset_draws \[`list()`\] of \[`data.table()`\]\cr
#'   Draws of directly modeled offset parameters for all ccmpp inputs as
#'   prepared in the `popReconstruct_*_draws()` functions.
#' @inheritParams extract_stan_draws
#' @param draw_col_name \[`character(1)`\]\cr
#'   Name of column in `offset_draws` \[`data.table()`\] that corresponds to
#'   draw index
#'
#' @return \[`list()`\] of \[`data.table()`\] containing draws of spline offsets
#'   parameters for all ccmpp inputs.
calculate_spline_offset_draws <- function(offset_draws,
                                          detailed_settings,
                                          draw_col_name = "draw") {

  spline_offset_draws <- lapply(names(offset_draws), function(comp) {
    offsets <- offset_draws[[comp]]
    spline_offsets <- offsets[, predict_spline_offset(
      dt = .SD,
      B_t = detailed_settings[[comp]][["B_t"]],
      B_a = detailed_settings[[comp]][["B_a"]],
      years = detailed_settings[[comp]][["years"]],
      ages = detailed_settings[[comp]][["ages"]]
    ), by = c(draw_col_name, if ("sex" %in% names(offsets)) "sex")]

    if ("year_start" %in% names(spline_offsets)) {
      hierarchyUtils::gen_end(
        dt = spline_offsets,
        id_cols = setdiff(names(offsets), c("value", "year_end", "age_end")),
        col_stem = "year",
        right_most_endpoint = max(offsets$year_end)
      )
    }

    if ("age_start" %in% names(spline_offsets)) {
      hierarchyUtils::gen_end(
        dt = spline_offsets,
        id_cols = setdiff(names(offsets), c("value", "year_end", "age_end")),
        col_stem = "age",
        right_most_endpoint = max(offsets$age_end)
      )
    }
    data.table::setkeyv(spline_offsets, setdiff(names(spline_offsets), "value"))
    return(spline_offsets)
  })

  names(spline_offset_draws) <- names(offset_draws)
  return(spline_offset_draws)
}

#' @title Helper function to predict spline offset for one ccmpp input component
#'
#' @param dt \[`data.table()`\]\cr
#'   Input offset estimates for the one ccmpp input component for one draw and
#'   sex.
#' @param B_t \[`matrix()`\]\cr
#'   B-spline basis matrix for the year dimension. By default is NULL but must
#'   be provided for all ccmpp inputs except 'baseline'.
#' @param B_a \[`matrix()`\]\cr
#'   B-spline basis matrix for the age dimension. By default is NULL but must
#'   be provided for all ccmpp inputs except 'srb'.
#' @param years \[`integer()`\]\cr
#'   The start of each calendar year interval for the ccmpp input component. By
#'   default is NULL but must be provided for all ccmpp inputs except 'baseline'.
#' @param ages \[`integer()`\]\cr
#'   The start of each age group for the ccmpp input component. By default is
#'   NULL but must be provided for all ccmpp inputs except 'srb'.
#'
#' @return \[`data.table()`\] containing draws of spline offsets for one ccmpp
#'   input.
predict_spline_offset <- function(dt, B_t = NULL, B_a = NULL,
                                  years = NULL, ages = NULL) {

  # reformat the offset from data.table to matrix with ages for rows and years
  # for columns
  form <- paste0(ifelse(is.null(B_a), ".", "age_start"), " ~ ",
                 ifelse(is.null(B_t), ".", "year_start"))
  offset_matrix <- dcast(dt, formula = stats::as.formula(form),
                         value.var = "value")
  offset_matrix[, ifelse(is.null(ages), ".", "age_start") := NULL]
  offset_matrix <- as.matrix(offset_matrix)

  # calculate spline offset matrix
  if (!is.null(B_a) & !is.null(B_t)) {
    spline_offset_matrix <- B_a %*% offset_matrix %*% t(B_t)
  } else if (!is.null(B_a)) {
    spline_offset_matrix <- B_a %*% offset_matrix
  } else if (!is.null(B_t)) {
    spline_offset_matrix <- offset_matrix %*% t(B_t)
  }

  if (!is.null(years)) colnames(spline_offset_matrix) <- years
  if (!is.null(ages)) {
    rownames(spline_offset_matrix) <- ages
  } else {
    rownames(spline_offset_matrix) <- 0
  }

  # reformat the spline offset matrix from matrix to data.table
  spline_offset <- demCore:::matrix_to_dt(
    mdt = spline_offset_matrix,
    gen_end_interval_col = FALSE,
    validate_arguments = FALSE
  )

  calendar_interval_input <- "year_end" %in% names(dt)
  if (!calendar_interval_input) {
    setnames(spline_offset, "year_start", "year")
  }

  return(spline_offset)
}

#' @title Helper function to predict ccmpp inputs estimates from spline offsets
#'   and initial estimates.
#'
#' @param spline_offset_draws \[`list()`\] of \[`data.table()`\]\cr
#'   Draws of spline offsets as returned by [`calculate_spline_offset_draws()`].
#' @inheritParams popReconstruct_fit
#' @inheritParams calculate_spline_offset_draws
#'
#' @return \[`list()`\] of \[`data.table()`\] containing draws of all ccmpp
#'   inputs.
calculate_ccmpp_input_draws <- function(spline_offset_draws,
                                        inputs,
                                        detailed_settings,
                                        value_col) {

  ccmpp_input_draws <- lapply(names(spline_offset_draws), function(comp) {

    spline_offset <- copy(spline_offset_draws[[comp]])
    setnames(spline_offset, "value", "spline_offset")

    input <- copy(inputs[[comp]])
    data.table::setnames(input, value_col, "initial")
    transform_dt(
      dt = input,
      value_col = "initial",
      transformation = detailed_settings[[comp]][["transformation"]],
      transformation_arguments = detailed_settings[[comp]][["transformation_arguments"]]
    )

    input_draws <- merge(
      x = spline_offset,
      y = input,
      all = TRUE,
      by = setdiff(names(input), "initial")
    )
    input_draws[, value := spline_offset + initial]
    transform_dt(
      dt = input_draws,
      value_col = "value",
      transformation = detailed_settings[[comp]][["inverse_transformation"]],
      transformation_arguments = detailed_settings[[comp]][["transformation_arguments"]]
    )

    input_draws[, c("spline_offset", "initial") := NULL]
    data.table::setkeyv(input_draws, setdiff(names(input_draws), "value"))
    return(input_draws)
  })
  names(ccmpp_input_draws) <- names(spline_offset_draws)
  return(ccmpp_input_draws)
}

#' @title Helper function to do ccmpp at the draw level using
#'   [`demCore::ccmpp()`]
#'
#' @param input_draws \[`list()`\] of \[`data.table()`\]\cr
#'   Draws of ccmpp inputs as returned by [`calculate_ccmpp_input_draws()`].
#' @inheritParams popReconstruct_fit
#' @param n_draws \[`integer()`\]\cr
#'   Number of draws that are included in each ccmpp input \[`data.table()`\].
#' @inheritParams  calculate_spline_offset_draws
#'
#' @return \[`data.table()`\] containing draws of projected population
#'   estimates.
ccmpp_draws <- function(input_draws,
                        settings,
                        n_draws,
                        draw_col_name) {

  pop_draws <- lapply(1:n_draws, function(i) {

    # get one draw of each of the input components
    input <- lapply(names(input_draws), function(comp) {
      one_draw_dt <- input_draws[[comp]][get(draw_col_name) == i]
      one_draw_dt[[draw_col_name]] <- NULL
      return(one_draw_dt)
    })
    names(input) <- names(input_draws)

    pop <- demCore::ccmpp(
      inputs = input,
      settings = settings,
      value_col = "value",
      assert_positive_pop = FALSE,
      validate_arguments = FALSE,
      gen_end_interval_col = FALSE
    )
    pop[, draw := i]
    return(pop)
  })
  pop_draws <- rbindlist(pop_draws)
  data.table::setnames(pop_draws, "draw", draw_col_name)
  data.table::setkeyv(pop_draws, setdiff(names(pop_draws), "value"))
  return(pop_draws)
}

#' @title Helper function for rbinding each \[`data.table()`\] nested within a
#'   \[`list()`\]
#'
#' @description Inputs multiple lists of data.tables that will be combined
#'   into one list of data.tables. data.tables are matched by the name of each
#'   element of each list. This is different than [`data.table::rbindlist()`]
#'   which inputs one list of data.tables and returns one combined data.table.
#'
#' @param ... multiple \[`list()`\] of many \[`data.table()`\]\cr
#'   as returned by [`popReconstruct_posterior_draws()`] or
#'   [`popReconstruct_prior_draws()`].
#'
#' @return one \[`list()`\] of many \[`data.table()`\]
#'
#' @family popReconstruct
#'
#' @examples
#' list_dt1 <- list(
#'   "population" = data.table::data.table(
#'     sex = c("female", "male"),
#'     value = 2,
#'     method = 1
#' ),
#'   "deaths" = data.table::data.table(
#'     sex = c("female", "male"),
#'     value = 4,
#'     method = 1
#'   )
#' )
#' list_dt2 <- list(
#'   "population" = data.table::data.table(
#'     sex = c("female", "male"),
#'     value = 2,
#'     method = 2
#'   ),
#'   "deaths" = data.table::data.table(
#'     sex = c("female", "male"),
#'     value = 4,
#'     method = 2
#'   )
#' )
#' combined_list_dt <- rbindlist_dts(list_dt1, list_dt2)
#'
#' @export
rbindlist_dts <- function(...) {

  dots <- list(...)

  # determine the named elements corresponding to the estimate tables
  measures <- lapply(1:length(dots), function(i) names(dots[[i]]))
  measures <- unique(unlist(measures))

  # combine together all estimates from different methods
  combined_results <- lapply(measures, function(measure){
    measure_results <- lapply(1:length(dots), function(i) {
      result <- copy(dots[[i]][[measure]])
      return(result)
    })
    measure_results <- data.table::rbindlist(measure_results, use.names = T, fill = T)

    # delete columns present in extracted stan results
    delete_columns <- c("chain", "chain_draw")
    if (any(delete_columns %in% names(measure_results))) measure_results[, c("chain", "chain_draw") := NULL]

    return(measure_results)
  })
  names(combined_results) <- measures
  return(combined_results)
}

#' @export
#' @rdname popReconstruct_fit
popReconstruct_summarize_draws <- function(draws,
                                           summarize_cols = c("chain", "chain_draw", "draw"),
                                           value_col = "value",
                                           ...) {

  # check `summarize_cols` argument
  assertive::assert_is_character(summarize_cols)

  # check `value_col` argument
  assertthat::assert_that(assertthat::is.string(value_col))

  # check `draws` argument
  assertthat::assert_that(
    class(draws) == "list",
    unique(sapply(draws, class)[1, ]) == "data.table",
    all(sapply(draws, function(dt) value_col %in% names(dt))),
    all(sapply(draws, function(dt) all(summarize_cols %in% names(dt)))),
    msg = paste0("`draws` must be a list of data.tables that each contains ",
                 "`summarize_cols` and `value_col`")
  )

  summaries <- lapply(draws, function(dt_draws) {
    id_cols <- setdiff(names(dt_draws), value_col)
    dt_summary <- demUtils::summarize_dt(
      dt = dt_draws,
      id_cols = id_cols,
      summarize_cols = summarize_cols,
      value_col = value_col,
      ...
    )
    return(dt_summary)
  })
  names(summaries) <- names(draws)
  return(summaries)
}
