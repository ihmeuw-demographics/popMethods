/// @file popReconstruct.hpp

#ifndef popReconstruct_hpp
#define popReconstruct_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// Probability density function of the inverse gamma distribution
template<class Type>
Type dinvGamma(Type x, Type shape, Type scale, int give_log=0) {
  if (give_log) {
    return shape * log(scale) - lgamma(asDouble(shape)) - (shape + Type(1)) * log(x) - (scale / x);
  } else {
    return pow(scale, shape) / tgamma(asDouble(shape)) * pow((Type(1) / x), (shape + Type(1))) * exp(-scale / x);
  }
}

template<class Type>
array<Type> bounded_inv_logit(array<Type> x, int domain_lower, int domain_upper) {

  array<Type> result(x.dim);
  result = exp(x) / (1 + exp(x));
  result = (result * (domain_upper - domain_lower)) + domain_lower;
  return result;
}

template<class Type>
matrix<Type> calculate_nSx(array<Type> mx, array<Type> non_terminal_ax,
                           array<Type> terminal_ax, int interval,
                           int A, int Y, int sexes) {

  array<Type> nSx(mx.dim);

  for (int s = 0; s < sexes; s++) {
    for (int y = 0; y < Y; y++) {
      vector<Type> qx(A + 1);
      vector<Type> px(A + 1);
      vector<Type> lx(A + 1);
      vector<Type> dx(A + 1);
      vector<Type> nLx(A + 1);
      vector<Type> Tx(A + 1);

      // calculate qx
      for (int a = 0; a < A; a++) {
        qx(a) = (interval * mx.col(s).col(y)(a)) / (1 + ((interval - non_terminal_ax.col(s).col(y)(a)) * mx.col(s).col(y)(a)));
      }
      qx(A + 1) = 1;

      // calculate px
      px = 1 - qx;

      // calculate lx
      lx(0) = 1;
      for (int a = 1; a < (A + 1); a++) {
        lx(a) = lx(a - 1) * px (a - 1);
      }

      // calculate dx
      for (int a = 0; a < A; a++) {
        dx(a) = lx(a) - lx(a + 1);
      }
      dx(A) = lx(A);

      // calculate nLx
      for (int a = 0; a < A; a++) {
        nLx(a) = (interval * lx(a + 1)) + (non_terminal_ax.col(s).col(y)(a) * dx(a));
      }
      nLx(A) = lx(A) / mx.col(s).col(y)(A);

      // calculate Tx
      for (int a1 = 0; a1 < (A + 1); a1++) {
        for (int a2 = a1; a2 < (A + 1); a2++) {
          Tx(a1) += nLx(a2);
        }
      }

      // calculate Sx
      nSx.col(s).col(y)(0) = nLx(0) / (interval * lx(0));
      for (int a = 1; a < A; a++) {
        nSx.col(s).col(y)(a) = nLx(a) / nLx(a - 1);
      }
      nSx.col(s).col(y)(A) = Tx(A) / Tx(A - 1);
    }
  }
  return(nSx);
}

template<class Type>
matrix<Type> make_leslie_matrix(Type srb, array<Type> asfr, array<Type> survival,
                                int interval, int A, int A_f, int A_f_offset, bool female) {

  // initialize leslie matrix to all zeroes
  matrix<Type> leslie;
  leslie.setZero(A, A);

  // first row includes asfr, srb, and birth survival (used to calculate youngest female population age group)
  if (female) {
    Type k = (1 / (1 + srb)) * survival(0) * 0.5 * interval;
    leslie(0, A_f_offset - 1) = k * asfr(0)* survival(A_f_offset); // fertility contribution from women aging into the youngest reproductive age group
    leslie(0, A_f_offset + A_f - 1) = k * asfr(A_f - 1); // fertility contribution from women starting in the oldest reproductive age group
    for (int a_f = 0; a_f < A_f - 1; a_f++) {
      int a = a_f + A_f_offset;
      leslie(0, a) = k * (asfr(a_f) + (asfr(a_f + 1) * survival(a + 1)));
    }
  }

  // other rows include survivorship ratios
  for (int a = 0; a < A - 1; a++) {
    leslie(a + 1, a) = survival(a + 1);
  }
  leslie(A - 1, A - 1) = survival(A);

  return leslie;
}

template<class Type>
array<Type> ccmpp(array<Type> srb, array<Type> asfr, array<Type> baseline,
                  array<Type> survival, array<Type> net_migration,
                  int sexes, int interval, int A, int A_f, int A_f_offset, int Y) {

  // initialize projected population counts
  array<Type> population(A, Y + 1, sexes);
  for (int s = 0; s < sexes; s++) {
    population.col(s).col(0) = baseline.col(s).col(0); // set first year equal to baseline
  }

  // project forward population one year at a time
  for (int y = 0; y < Y; y++) {
    // project female population forward one projection period
    matrix<Type> half_net_migration_count_f = population.col(0).col(y) * net_migration.col(0).col(y) * 0.5;
    matrix<Type> leslie_f = make_leslie_matrix(srb(0, y), asfr.col(y), survival.col(0).col(y), interval, A, A_f, A_f_offset, true);
    population.col(0).col(y + 1) = (leslie_f * (matrix<Type>(population.col(0).col(y)) + half_net_migration_count_f) + half_net_migration_count_f);

    // project male population forward one projection period
    if (sexes == 2) {
      matrix<Type> half_net_migration_count_m = population.col(1).col(y) * net_migration.col(1).col(y) * 0.5;
      matrix<Type> leslie_m = make_leslie_matrix(srb(0, y), asfr.col(y), survival.col(1).col(y), interval, A, A_f, A_f_offset, false);
      population.col(1).col(y + 1) = (leslie_m * (matrix<Type>(population.col(1).col(y)) + half_net_migration_count_m) + half_net_migration_count_m);

      // back-calculate total births in projection period
      Type total_births = population(0, y + 1, 0) / (survival(0, y, 0) * (1 / (1 + srb(0, y))));
      // calculate youngest male population age group
      Type young_male_pop = total_births * survival(0, y, 1) * (srb(0, y) / (1 + srb(0, y)));
      // add onto migrants that are already calculated
      population(0, y + 1, 1) = population(0, y + 1, 1) + young_male_pop;
    }
  }

  return(population);
}

template<class Type>
vector<Type> aggregate(array<Type> population, array<int> input_pop_year_index,
                       array<int> input_pop_age_index, array<int> input_pop_sex_index,
                       int N_pop) {

  vector<Type> population_aggregated(N_pop);
  for (int i = 0; i < N_pop; i++) {
    Type sum_pop = 0;
    for (int y = input_pop_year_index(i, 0); y < input_pop_year_index(i, 1); y++) {
      for (int a = input_pop_age_index(i, 0); a < input_pop_age_index(i, 1); a++) {
        for (int s = input_pop_sex_index(i, 0); s < input_pop_sex_index(i, 1); s++) {
          sum_pop += population(a, y, s);
        }
      }
    }
    population_aggregated(i) = sum_pop;
  }

  return(population_aggregated);
}

/// Negative log-likelihood of the popReconstruct model.
template<class Type>
Type popReconstruct(objective_function<Type>* obj) {

  ////// data section

  // general setup variables
  DATA_INTEGER(sexes); // number of sexes (1 for female only, 2 for female and male)
  DATA_INTEGER(interval); // age and year interval
  DATA_INTEGER(A); // number of age groups to estimate for
  DATA_INTEGER(A_f); // number of reproductive age groups to estimate for
  DATA_INTEGER(A_f_offset); // number of younger age groups that are not included in the reproductive ages
  DATA_INTEGER(Y); // number of year intervals in projection period
  DATA_INTEGER(Y_p); // number of population data years (not including the baseline year)
  DATA_INTEGER(estimate_survival); // whether to estimate survival (or mx and ax separately)
  DATA_INTEGER(estimate_net_migration); // whether to estimate net migration (or immigration and emigration separately)

  // inverse gamma alpha and beta hyperparameters for ccmpp inputs
  DATA_SCALAR(alpha_srb);
  DATA_SCALAR(beta_srb);
  DATA_SCALAR(alpha_asfr);
  DATA_SCALAR(beta_asfr);
  DATA_SCALAR(alpha_population);
  DATA_SCALAR(beta_population);
  DATA_SCALAR(alpha_survival);
  DATA_SCALAR(beta_survival);
  DATA_SCALAR(alpha_mx);
  DATA_SCALAR(beta_mx);
  DATA_SCALAR(alpha_non_terminal_ax);
  DATA_SCALAR(beta_non_terminal_ax);
  DATA_SCALAR(alpha_terminal_ax);
  DATA_SCALAR(beta_terminal_ax);
  DATA_SCALAR(alpha_net_migration);
  DATA_SCALAR(beta_net_migration);
  DATA_SCALAR(alpha_immigration);
  DATA_SCALAR(beta_immigration);
  DATA_SCALAR(alpha_emigration);
  DATA_SCALAR(beta_emigration);

  // ccmpp input initial estimates
  DATA_ARRAY(input_log_srb); // input sex ratio at birth
  DATA_ARRAY(input_log_asfr); // input age-specific fertility rates
  DATA_ARRAY(input_log_baseline); // input baseline populations
  DATA_ARRAY(input_logit_survival); // input survival proportions
  DATA_ARRAY(input_log_mx); // input mortality rate
  DATA_ARRAY(input_bounded_logit_non_terminal_ax); // input ax for non-terminal age groups
  DATA_ARRAY(input_log_terminal_ax); // input ax for terminal age group
  DATA_ARRAY(input_net_migration); // input net migration proportions
  DATA_ARRAY(input_log_immigration); // input immigration proportions
  DATA_ARRAY(input_log_emigration); // input emigration proportions

  // B-spline linear basis functions
  DATA_INTEGER(N_k_t_srb);
  DATA_MATRIX(B_t_srb);

  DATA_INTEGER(N_k_t_asfr);
  DATA_MATRIX(B_t_asfr);
  DATA_INTEGER(N_k_a_asfr);
  DATA_MATRIX(B_a_asfr);

  DATA_INTEGER(N_k_a_baseline);
  DATA_MATRIX(B_a_baseline);

  DATA_INTEGER(N_k_t_survival);
  DATA_MATRIX(B_t_survival);
  DATA_INTEGER(N_k_a_survival);
  DATA_MATRIX(B_a_survival);

  DATA_INTEGER(N_k_t_mx);
  DATA_MATRIX(B_t_mx);
  DATA_INTEGER(N_k_a_mx);
  DATA_MATRIX(B_a_mx);

  DATA_INTEGER(N_k_t_non_terminal_ax);
  DATA_MATRIX(B_t_non_terminal_ax);
  DATA_INTEGER(N_k_a_non_terminal_ax);
  DATA_MATRIX(B_a_non_terminal_ax);

  DATA_INTEGER(N_k_t_terminal_ax);
  DATA_MATRIX(B_t_terminal_ax);

  DATA_INTEGER(N_k_t_net_migration);
  DATA_MATRIX(B_t_net_migration);
  DATA_INTEGER(N_k_a_net_migration);
  DATA_MATRIX(B_a_net_migration);

  DATA_INTEGER(N_k_t_immigration);
  DATA_MATRIX(B_t_immigration);
  DATA_INTEGER(N_k_a_immigration);
  DATA_MATRIX(B_a_immigration);

  DATA_INTEGER(N_k_t_emigration);
  DATA_MATRIX(B_t_emigration);
  DATA_INTEGER(N_k_a_emigration);
  DATA_MATRIX(B_a_emigration);

  // population data
  // population data
  DATA_INTEGER(N_pop); // number of year-age-sex specific population data points
  DATA_VECTOR(input_log_pop_value); // log population value for each year-age-sex specific data point
  DATA_VECTOR(input_pop_weight); // population variance weighting value for each year-age-sex specific data point
  DATA_IARRAY(input_pop_year_index); // start (inclusive) and end (exclusive) indices for year in full 'population' object below
  DATA_IARRAY(input_pop_age_index); // start (inclusive) and end (exclusive) indices for age in full 'population' object below
  DATA_IARRAY(input_pop_sex_index); // start (inclusive) and end (exclusive) indices for sex in full 'population' object below

  ////// parameter section

  // log variance parameters
  PARAMETER(log_sigma2_srb);
  PARAMETER(log_sigma2_asfr);
  PARAMETER(log_sigma2_population);
  PARAMETER(log_sigma2_survival);
  PARAMETER(log_sigma2_mx);
  PARAMETER(log_sigma2_non_terminal_ax);
  PARAMETER(log_sigma2_terminal_ax);
  PARAMETER(log_sigma2_net_migration);
  PARAMETER(log_sigma2_immigration);
  PARAMETER(log_sigma2_emigration);

  // offset parameters
  PARAMETER_ARRAY(offset_log_srb);
  PARAMETER_ARRAY(offset_log_asfr);
  PARAMETER_ARRAY(offset_log_baseline);
  PARAMETER_ARRAY(offset_logit_survival);
  PARAMETER_ARRAY(offset_log_mx);
  PARAMETER_ARRAY(offset_bounded_logit_non_terminal_ax);
  PARAMETER_ARRAY(offset_log_terminal_ax);
  PARAMETER_ARRAY(offset_net_migration);
  PARAMETER_ARRAY(offset_log_immigration);
  PARAMETER_ARRAY(offset_log_emigration);


  ////// transformed parameters section

  // variance parameters
  Type sigma2_srb = exp(log_sigma2_srb);
  Type sigma2_asfr = exp(log_sigma2_asfr);
  Type sigma2_population = exp(log_sigma2_population);
  Type sigma2_survival = exp(log_sigma2_survival);
  Type sigma2_mx = exp(log_sigma2_mx);
  Type sigma2_non_terminal_ax = exp(log_sigma2_non_terminal_ax);
  Type sigma2_terminal_ax = exp(log_sigma2_terminal_ax);
  Type sigma2_net_migration = exp(log_sigma2_net_migration);
  Type sigma2_immigration = exp(log_sigma2_immigration);
  Type sigma2_emigration = exp(log_sigma2_emigration);

  // spline offset parameters
  array<Type> spline_offset_log_srb(input_log_srb.dim);
  array<Type> spline_offset_log_asfr(input_log_asfr.dim);
  array<Type> spline_offset_log_baseline(input_log_baseline.dim);
  array<Type> spline_offset_logit_survival(input_logit_survival.dim);
  array<Type> spline_offset_log_mx(input_log_mx.dim);
  array<Type> spline_offset_bounded_logit_non_terminal_ax(input_bounded_logit_non_terminal_ax.dim);
  array<Type> spline_offset_log_terminal_ax(input_log_terminal_ax.dim);
  array<Type> spline_offset_net_migration(input_net_migration.dim);
  array<Type> spline_offset_log_immigration(input_log_immigration.dim);
  array<Type> spline_offset_log_emigration(input_log_emigration.dim);

  // untransformed ccmp input parameters
  array<Type> srb(input_log_srb.dim);
  array<Type> asfr(input_log_asfr.dim);
  array<Type> baseline(input_log_baseline.dim);
  array<Type> survival(input_logit_survival.dim);
  array<Type> mx(input_log_mx.dim);
  array<Type> non_terminal_ax(input_bounded_logit_non_terminal_ax.dim);
  array<Type> terminal_ax(input_log_terminal_ax.dim);
  array<Type> net_migration(input_net_migration.dim);
  array<Type> immigration(input_log_immigration.dim);
  array<Type> emigration(input_log_emigration.dim);

  // calculate spline offsets
  spline_offset_log_srb = matrix<Type>(offset_log_srb) * B_t_srb.transpose();
  spline_offset_log_asfr = B_a_asfr * matrix<Type>(offset_log_asfr) * B_t_asfr.transpose();
  // need to do sex specific offsets like this since TMB collapses array dimensions
  for (int s = 0; s < sexes; s++) {
    matrix<Type> sex_specific_offset_log_baseline_matrix(N_k_a_baseline, 1);
    matrix<Type> sex_specific_offset_logit_survival_matrix(N_k_a_survival, N_k_t_survival);
    matrix<Type> sex_specific_offset_log_mx_matrix(N_k_a_mx, N_k_t_mx);
    matrix<Type> sex_specific_offset_bounded_logit_non_terminal_ax_matrix(N_k_a_non_terminal_ax, N_k_t_non_terminal_ax);
    matrix<Type> sex_specific_offset_log_terminal_ax_matrix(1, N_k_t_terminal_ax);
    matrix<Type> sex_specific_offset_net_migration_matrix(N_k_a_net_migration, N_k_t_net_migration);
    matrix<Type> sex_specific_offset_log_immigration_matrix(N_k_a_immigration, N_k_t_immigration);
    matrix<Type> sex_specific_offset_log_emigration_matrix(N_k_a_emigration, N_k_t_emigration);

    for (int a = 0; a < N_k_a_baseline; a++) {
      sex_specific_offset_log_baseline_matrix(a, 0) = offset_log_baseline(a, 0, s);
    }
    spline_offset_log_baseline.col(s) = B_a_baseline * sex_specific_offset_log_baseline_matrix;

    // calculate either the survival spline offsets or mx & ax spline offsets
    if (estimate_survival) {
      for (int y = 0; y < N_k_t_survival; y++) {
        for (int a = 0; a < N_k_a_survival; a++) {
          sex_specific_offset_logit_survival_matrix(a, y) = offset_logit_survival(a, y, s);
        }
      }
      spline_offset_logit_survival.col(s) = B_a_survival * sex_specific_offset_logit_survival_matrix * B_t_survival.transpose();
    } else {
      for (int y = 0; y < N_k_t_mx; y++) {
        for (int a = 0; a < N_k_a_mx; a++) {
          sex_specific_offset_log_mx_matrix(a, y) = offset_log_mx(a, y, s);
        }
      }
      spline_offset_log_mx.col(s) = B_a_mx * sex_specific_offset_log_mx_matrix * B_t_mx.transpose();

      for (int y = 0; y < N_k_t_non_terminal_ax; y++) {
        for (int a = 0; a < N_k_a_non_terminal_ax; a++) {
          sex_specific_offset_bounded_logit_non_terminal_ax_matrix(a, y) = offset_bounded_logit_non_terminal_ax(a, y, s);
        }
      }
      spline_offset_bounded_logit_non_terminal_ax.col(s) = B_a_non_terminal_ax * sex_specific_offset_bounded_logit_non_terminal_ax_matrix * B_t_non_terminal_ax.transpose();

      for (int y = 0; y < N_k_t_terminal_ax; y++) {
        sex_specific_offset_log_terminal_ax_matrix(0, y) = offset_log_terminal_ax(0, y, s);
      }
      spline_offset_log_terminal_ax.col(s) = sex_specific_offset_log_terminal_ax_matrix * B_t_terminal_ax.transpose();
    }

    // calculate either the net migration spline offsets or immigration & emigration spline offsets
    if (estimate_net_migration) {
      for (int y = 0; y < N_k_t_net_migration; y++) {
        for (int a = 0; a < N_k_a_net_migration; a++) {
          sex_specific_offset_net_migration_matrix(a, y) = offset_net_migration(a, y, s);
        }
      }
      spline_offset_net_migration.col(s) = B_a_net_migration * sex_specific_offset_net_migration_matrix * B_t_net_migration.transpose();
    } else {
      for (int y = 0; y < N_k_t_immigration; y++) {
        for (int a = 0; a < N_k_a_immigration; a++) {
          sex_specific_offset_log_immigration_matrix(a, y) = offset_log_immigration(a, y, s);
        }
      }
      spline_offset_log_immigration.col(s) = B_a_immigration * sex_specific_offset_log_immigration_matrix * B_t_immigration.transpose();

      for (int y = 0; y < N_k_t_emigration; y++) {
        for (int a = 0; a < N_k_a_emigration; a++) {
          sex_specific_offset_log_emigration_matrix(a, y) = offset_log_emigration(a, y, s);
        }
      }
      spline_offset_log_emigration.col(s) = B_a_emigration * sex_specific_offset_log_emigration_matrix * B_t_emigration.transpose();
    }
  }

  // calculate untransformed ccmpp input parameters
  srb = exp(input_log_srb + spline_offset_log_srb);
  asfr = exp(input_log_asfr + spline_offset_log_asfr);
  baseline = exp(input_log_baseline + spline_offset_log_baseline);
  if (estimate_survival) {
    survival = bounded_inv_logit(input_logit_survival + spline_offset_logit_survival, 0, 1);
  } else {
    mx = exp(input_log_mx + spline_offset_log_mx);
    non_terminal_ax = bounded_inv_logit(input_bounded_logit_non_terminal_ax + spline_offset_bounded_logit_non_terminal_ax, 0, interval);
    terminal_ax = exp(input_log_terminal_ax + spline_offset_log_terminal_ax);
    survival = calculate_nSx(mx, non_terminal_ax, terminal_ax, interval, A, Y, sexes);
  }
  if (estimate_net_migration) {
    net_migration = input_net_migration + spline_offset_net_migration;
  } else {
    immigration = exp(input_log_immigration + spline_offset_log_immigration);
    emigration = exp(input_log_emigration + spline_offset_log_emigration);
    net_migration = immigration - emigration;
  }

  // Level 2 (ccmpp)
  array<Type> population(A, Y + 1, sexes);
  // project population with most detailed age groups
  population = ccmpp(srb, asfr, baseline, survival, net_migration,
                     sexes, interval, A, A_f, A_f_offset, Y);
  // aggregate to census age groups
  vector<Type> population_input_groups(N_pop);
  population_input_groups = aggregate(population, input_pop_year_index,
                                      input_pop_age_index, input_pop_sex_index,
                                      N_pop);


  ////// model section

  Type nll = 0;

  // LEVEL 4 (informative prior distributions for variance parameters)
  nll -= dinvGamma(sigma2_srb, alpha_srb, beta_srb, true);
  nll -= dinvGamma(sigma2_asfr, alpha_asfr, beta_asfr, true);
  nll -= dinvGamma(sigma2_population, alpha_population, beta_population, true);
  nll -= dinvGamma(sigma2_survival, alpha_survival, beta_survival, true);
  nll -= dinvGamma(sigma2_mx, alpha_mx, beta_mx, true);
  nll -= dinvGamma(sigma2_non_terminal_ax, alpha_non_terminal_ax, beta_non_terminal_ax, true);
  nll -= dinvGamma(sigma2_terminal_ax, alpha_terminal_ax, beta_terminal_ax, true);
  nll -= dinvGamma(sigma2_net_migration, alpha_net_migration, beta_net_migration, true);
  nll -= dinvGamma(sigma2_immigration, alpha_immigration, beta_immigration, true);
  nll -= dinvGamma(sigma2_emigration, alpha_emigration, beta_emigration, true);

  // Jacobian adjustment for variances (this is automatically handled by Stan)
  // See https://www.rpubs.com/kaz_yos/stan_jacobian for example
  // See https://github.com/colemonnahan/hmc_tests/blob/master/models/swallows/swallows.cpp for example
  nll -= log_sigma2_srb + log_sigma2_asfr + log_sigma2_population + log_sigma2_survival + log_sigma2_mx + log_sigma2_non_terminal_ax + log_sigma2_terminal_ax + log_sigma2_net_migration + log_sigma2_immigration + log_sigma2_emigration;

  // LEVEL 3 (model initial estimates ccmpp of ccmpp inputs)
  nll -= dnorm(vector<Type>(offset_log_srb), 0, sqrt(sigma2_srb), true).sum();
  nll -= dnorm(vector<Type>(offset_log_asfr), 0, sqrt(sigma2_asfr), true).sum();
  nll -= dnorm(vector<Type>(offset_log_baseline), 0, sqrt(sigma2_population), true).sum();
  nll -= dnorm(vector<Type>(offset_logit_survival), 0, sqrt(sigma2_survival), true).sum();
  nll -= dnorm(vector<Type>(offset_log_mx), 0, sqrt(sigma2_mx), true).sum();
  nll -= dnorm(vector<Type>(offset_bounded_logit_non_terminal_ax), 0, sqrt(sigma2_non_terminal_ax), true).sum();
  nll -= dnorm(vector<Type>(offset_log_terminal_ax), 0, sqrt(sigma2_terminal_ax), true).sum();
  nll -= dnorm(vector<Type>(offset_net_migration), 0, sigma2_net_migration, true).sum();
  nll -= dnorm(vector<Type>(offset_log_immigration), 0, sigma2_immigration, true).sum();
  nll -= dnorm(vector<Type>(offset_log_emigration), 0, sigma2_emigration, true).sum();

  // LEVEL 1 (model the census counts)
  nll -= dnorm(input_log_pop_value, log(population_input_groups), sqrt(sigma2_population) * input_pop_weight, true).sum();

  return nll;
}


#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
