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
  DATA_IVECTOR(pop_data_years_index); // columns in population object that correspond to population data years

  // inverse gamma alpha and beta hyperparameters for ccmpp inputs
  DATA_SCALAR(alpha_srb);
  DATA_SCALAR(beta_srb);
  DATA_SCALAR(alpha_asfr);
  DATA_SCALAR(beta_asfr);
  DATA_SCALAR(alpha_population);
  DATA_SCALAR(beta_population);
  DATA_SCALAR(alpha_survival);
  DATA_SCALAR(beta_survival);
  DATA_SCALAR(alpha_net_migration);
  DATA_SCALAR(beta_net_migration);

  // ccmpp input initial estimates
  DATA_ARRAY(input_log_srb); // input sex ratio at birth
  DATA_ARRAY(input_log_asfr); // input age-specific fertility rates
  DATA_ARRAY(input_log_baseline); // input baseline populations
  DATA_ARRAY(input_logit_survival); // input survival proportions
  DATA_ARRAY(input_net_migration); // input net migration proportions

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

  DATA_INTEGER(N_k_t_net_migration);
  DATA_MATRIX(B_t_net_migration);
  DATA_INTEGER(N_k_a_net_migration);
  DATA_MATRIX(B_a_net_migration);

  // population data
  DATA_ARRAY(input_log_population_data); // population data


  ////// parameter section

  // log variance parameters
  PARAMETER(log_sigma2_srb);
  PARAMETER(log_sigma2_asfr);
  PARAMETER(log_sigma2_population);
  PARAMETER(log_sigma2_survival);
  PARAMETER(log_sigma2_net_migration);

  // offset parameters
  PARAMETER_ARRAY(offset_log_srb);
  PARAMETER_ARRAY(offset_log_asfr);
  PARAMETER_ARRAY(offset_log_baseline);
  PARAMETER_ARRAY(offset_logit_survival);
  PARAMETER_ARRAY(offset_net_migration);


  ////// transformed parameters section

  // variance parameters
  Type sigma2_srb = exp(log_sigma2_srb);
  Type sigma2_asfr = exp(log_sigma2_asfr);
  Type sigma2_population = exp(log_sigma2_population);
  Type sigma2_survival = exp(log_sigma2_survival);
  Type sigma2_net_migration = exp(log_sigma2_net_migration);

  // spline offset parameters
  array<Type> spline_offset_log_srb(input_log_srb.dim);
  array<Type> spline_offset_log_asfr(input_log_asfr.dim);
  array<Type> spline_offset_log_baseline(input_log_baseline.dim);
  array<Type> spline_offset_logit_survival(input_logit_survival.dim);
  array<Type> spline_offset_net_migration(input_net_migration.dim);

  // untransformed ccmp input parameters
  array<Type> srb(input_log_srb.dim);
  array<Type> asfr(input_log_asfr.dim);
  array<Type> baseline(input_log_baseline.dim);
  array<Type> survival(input_logit_survival.dim);
  array<Type> net_migration(input_net_migration.dim);

  // calculate spline offsets
  spline_offset_log_srb = matrix<Type>(offset_log_srb) * B_t_srb.transpose();
  spline_offset_log_asfr = B_a_asfr * matrix<Type>(offset_log_asfr) * B_t_asfr.transpose();
  // need to do sex specific offsets like this since TMB collapses array dimensions
  for (int s = 0; s < sexes; s++) {
    matrix<Type> sex_specific_offset_log_baseline_matrix(N_k_a_baseline, 1);
    matrix<Type> sex_specific_offset_logit_survival_matrix(N_k_a_survival, N_k_t_survival);
    matrix<Type> sex_specific_offset_net_migration_matrix(N_k_a_net_migration, N_k_t_net_migration);

    for (int a = 0; a < N_k_a_baseline; a++) {
      sex_specific_offset_log_baseline_matrix(a, 0) = offset_log_baseline(a, 0, s);
    }
    spline_offset_log_baseline.col(s) = B_a_baseline * sex_specific_offset_log_baseline_matrix;

    for (int y = 0; y < N_k_t_survival; y++) {
      for (int a = 0; a < N_k_a_survival; a++) {
        sex_specific_offset_logit_survival_matrix(a, y) = offset_logit_survival(a, y, s);
      }
    }
    spline_offset_logit_survival.col(s) = B_a_survival * sex_specific_offset_logit_survival_matrix * B_t_survival.transpose();

    for (int y = 0; y < N_k_t_net_migration; y++) {
      for (int a = 0; a < N_k_a_net_migration; a++) {
        sex_specific_offset_net_migration_matrix(a, y) = offset_net_migration(a, y, s);
      }
    }
    spline_offset_net_migration.col(s) = B_a_net_migration * sex_specific_offset_net_migration_matrix * B_t_net_migration.transpose();
  }

  // calculate untransformed ccmpp input parameters
  srb = exp(input_log_srb + spline_offset_log_srb);
  asfr = exp(input_log_asfr + spline_offset_log_asfr);
  baseline = exp(input_log_baseline + spline_offset_log_baseline);
  survival = 1 / (1 + exp(-1 * (input_logit_survival + spline_offset_logit_survival)));
  net_migration = input_net_migration + spline_offset_net_migration;

  // Level 2 (ccmpp)
  array<Type> population(A, Y + 1, sexes);
  population = ccmpp(srb, asfr, baseline, survival, net_migration,
                     sexes, interval, A, A_f, A_f_offset, Y);


  ////// model section

  Type nll = 0;

  // LEVEL 4 (informative prior distributions for variance parameters)
  nll -= dinvGamma(sigma2_srb, alpha_srb, beta_srb, true);
  nll -= dinvGamma(sigma2_asfr, alpha_asfr, beta_asfr, true);
  nll -= dinvGamma(sigma2_population, alpha_population, beta_population, true);
  nll -= dinvGamma(sigma2_survival, alpha_survival, beta_survival, true);
  nll -= dinvGamma(sigma2_net_migration, alpha_net_migration, beta_net_migration, true);

  // Jacobian adjustment for variances (this is automatically handled by Stan)
  // See https://www.rpubs.com/kaz_yos/stan_jacobian for example
  // See https://github.com/colemonnahan/hmc_tests/blob/master/models/swallows/swallows.cpp for example
  nll -= log_sigma2_srb + log_sigma2_asfr + log_sigma2_population + log_sigma2_survival + log_sigma2_net_migration;

  // LEVEL 3 (model initial estimates ccmpp of ccmpp inputs)
  nll -= dnorm(vector<Type>(offset_log_srb), 0, sqrt(sigma2_srb), true).sum();
  nll -= dnorm(vector<Type>(offset_log_asfr), 0, sqrt(sigma2_asfr), true).sum();
  nll -= dnorm(vector<Type>(offset_log_baseline), 0, sqrt(sigma2_population), true).sum();
  nll -= dnorm(vector<Type>(offset_logit_survival), 0, sqrt(sigma2_survival), true).sum();
  nll -= dnorm(vector<Type>(offset_net_migration), 0, sigma2_net_migration, true).sum();

  // LEVEL 1 (model the census counts)
  for (int s = 0; s < sexes; s++) {
    for (int y_p = 0; y_p < Y_p; y_p++) {
      int y = pop_data_years_index(y_p);
      for (int a = 0; a < A; a++) {
        nll -= dnorm(input_log_population_data(a, y_p, s), log(population(a, y, s)), sigma2_population, true);
      }
    }
  }

  return nll;
}


#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
