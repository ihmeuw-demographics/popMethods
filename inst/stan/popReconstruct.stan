functions {
  matrix make_leslie_matrix(real srb, vector asfr, vector survival, int interval,
                            int A, int A_f, int A_f_offset, int female) {

    matrix[A, A] leslie = rep_matrix(0, A, A); // initialize leslie matrix to all zeroes

    // first row includes asfr, srb, and birth survival (used to calculate youngest female population age group)
    if (female) {
      real k = (1 / (1 + srb)) * survival[1] * 0.5 * interval;
      leslie[1, A_f_offset] = k * asfr[1] * survival[A_f_offset + 1]; // fertility contribution from women aging into the youngest reproductive age group
      leslie[1, A_f_offset + A_f] = k * asfr[A_f]; // fertility contribution from women starting in the oldest reproductive age group
      for (a_f in 1:(A_f - 1)) {
        int a = a_f + A_f_offset;
        leslie[1, a] = k * (asfr[a_f] + (asfr[a_f + 1] * survival[a + 1]));
      }
    }

    // other rows include survivorship ratios
    for (a in 1:(A - 1)) {
      leslie[a + 1, a] = survival[a + 1];
    }
    leslie[A, A] = survival[A + 1];

    return leslie;
  }

  matrix[] ccmpp(matrix srb, matrix asfr, matrix[] baseline, matrix[] survival,
                 matrix[] net_migration, int sexes, int interval,
                 int A, int A_f, int A_f_offset, int Y) {

    // initialize projected population counts
    matrix[A, Y + 1] population[sexes];
    for (s in 1:sexes) {
      population[s, , 1] = baseline[s, , 1]; // set first year equal to baseline
    }

    // project forward population one year at a time
    for (y in 1:Y) {
      // project female population forward one projection period
      vector[A] half_net_migration_count_f = population[1, , y] .* net_migration[1, , y] * 0.5;
      matrix[A, A] leslie_f = make_leslie_matrix(srb[1, y], asfr[, y], survival[1, , y], interval, A, A_f, A_f_offset, 1);
      population[1, , y + 1] = (leslie_f * (population[1, , y] + half_net_migration_count_f)) + half_net_migration_count_f;

      // project male population forward one projection period
      if (sexes == 2) {
        real total_births;
        real young_male_pop;

        vector[A] half_net_migration_count_m = population[2, , y] .* net_migration[2, , y] * 0.5;
        matrix[A, A] leslie_m = make_leslie_matrix(srb[1, y], asfr[, y], survival[2, , y], interval, A, A_f, A_f_offset, 0);
        population[2, , y + 1] = (leslie_m * (population[2, , y] + half_net_migration_count_m)) + half_net_migration_count_m;

        // back-calculate total births in projection period
        total_births = population[1, 1, y + 1] / (survival[1, 1, y] * (1 / (1 + srb[1, y])));
        // calculate youngest male population age group
        young_male_pop = total_births * survival[2, 1, y] * (srb[1, y] / (1 + srb[1, y]));
        // add onto migrants that are already calculated
        population[2, 1, y + 1] = population[2, 1, y + 1] + young_male_pop;
      }
    }

    return(population);
  }

  real[] aggregate(matrix[] population, int[, ] input_pop_year_index,
                   int[, ] input_pop_age_index, int[, ] input_pop_sex_index,
                   int N_pop) {

    real population_aggregated[N_pop];
    for (i in 1:N_pop) {
      real sum_pop = 0;
      for (y in input_pop_year_index[i, 1]:(input_pop_year_index[i, 2] - 1)) {
        for (a in input_pop_age_index[i, 1]:(input_pop_age_index[i, 2] - 1)) {
          for (s in input_pop_sex_index[i, 1]:(input_pop_sex_index[i, 2] - 1)) {
            sum_pop += population[s, a, y];
          }
        }
      }
      population_aggregated[i] = sum_pop;
    }

    return population_aggregated;
  }
}
data {
  // general setup variables
  int<lower = 1, upper = 2> sexes; // number of sexes (1 for female only, 2 for female and male)
  int interval; // age and year interval
  int A; // number of age groups to estimate for
  int A_f; // number of reproductive age groups to estimate for
  int A_f_offset; // number of younger age groups that are not included in the reproductive ages
  int Y; // number of year intervals in projection period
  int Y_p; // number of population data years (not including the baseline year)

  // whether to estimate or fix certain model components
  int<lower = 0, upper = 1> estimate_srb;
  int<lower = 0, upper = 1> estimate_asfr;
  int<lower = 0, upper = 1> estimate_baseline;
  int<lower = 0, upper = 1> estimate_survival;
  int<lower = 0, upper = 1> estimate_net_migration;
  int<lower = 0, upper = 1> estimate_immigration;
  int<lower = 0, upper = 1> estimate_emigration;

  // inverse gamma alpha and beta hyperparameters for ccmpp inputs
  real alpha_srb;
  real beta_srb;
  real alpha_asfr;
  real beta_asfr;
  real alpha_population;
  real beta_population;
  real alpha_survival;
  real beta_survival;
  real alpha_net_migration;
  real beta_net_migration;
  real alpha_immigration;
  real beta_immigration;
  real alpha_emigration;
  real beta_emigration;

  // ccmpp input initial estimates
  matrix[1, Y] input_log_srb; // input sex ratio at birth
  matrix[A_f, Y] input_log_asfr; // input age-specific fertility rates
  matrix[A, 1] input_log_baseline[sexes]; // input baseline populations
  matrix[A + 1, Y] input_logit_survival[sexes]; // input survival proportions
  matrix[A, Y] input_net_migration[sexes]; // input net migration proportions
  matrix[A, Y] input_log_immigration[sexes]; // input immigration proportions
  matrix[A, Y] input_log_emigration[sexes]; // input emigration proportions

  // B-spline linear basis functions
  int N_k_t_srb;
  matrix[Y, N_k_t_srb] B_t_srb;

  int N_k_t_asfr;
  matrix[Y, N_k_t_asfr] B_t_asfr;
  int N_k_a_asfr;
  matrix[A_f, N_k_a_asfr] B_a_asfr;

  int N_k_a_baseline;
  matrix[A, N_k_a_baseline] B_a_baseline;

  int N_k_t_survival;
  matrix[Y, N_k_t_survival] B_t_survival;
  int N_k_a_survival;
  matrix[A + 1, N_k_a_survival] B_a_survival;

  int N_k_t_net_migration;
  matrix[Y, N_k_t_net_migration] B_t_net_migration;
  int N_k_a_net_migration;
  matrix[A, N_k_a_net_migration] B_a_net_migration;

  int N_k_t_immigration;
  matrix[Y, N_k_t_immigration] B_t_immigration;
  int N_k_a_immigration;
  matrix[A, N_k_a_immigration] B_a_immigration;

  int N_k_t_emigration;
  matrix[Y, N_k_t_emigration] B_t_emigration;
  int N_k_a_emigration;
  matrix[A, N_k_a_emigration] B_a_emigration;

  // population data
  int N_pop; // number of year-age-sex specific population data points
  real input_log_pop_value[N_pop]; // log population value for each year-age-sex specific data point
  vector[N_pop] input_pop_weight; // population variance weighting value for each year-age-sex specific data point
  int<lower = 1, upper = Y + 2>  input_pop_year_index[N_pop, 2]; // start (inclusive) and end (exclusive) indices for year in full 'population' object below
  int<lower = 1, upper = A + 1> input_pop_age_index[N_pop, 2]; // start (inclusive) and end (exclusive) indices for age in full 'population' object below
  int<lower = 1, upper = sexes + 1> input_pop_sex_index[N_pop, 2]; // start (inclusive) and end (exclusive) indices for sex in full 'population' object below
}
parameters {
  // variance parameters
  real<lower = 0> sigma2_srb[estimate_srb];
  real<lower = 0> sigma2_asfr[estimate_asfr];
  real<lower = 0> sigma2_population;
  real<lower = 0> sigma2_survival[estimate_survival];
  real<lower = 0> sigma2_net_migration[estimate_net_migration];
  real<lower = 0> sigma2_immigration[estimate_immigration];
  real<lower = 0> sigma2_emigration[estimate_emigration];

  // offset parameters
  matrix[1, N_k_t_srb] offset_log_srb[estimate_srb];
  matrix[N_k_a_asfr, N_k_t_asfr] offset_log_asfr[estimate_asfr];
  matrix[N_k_a_baseline, 1] offset_log_baseline[estimate_baseline, sexes];
  matrix[N_k_a_survival, N_k_t_survival] offset_logit_survival[estimate_survival, sexes];
  matrix[N_k_a_net_migration, N_k_t_net_migration] offset_net_migration[estimate_net_migration, sexes];
  matrix[N_k_a_immigration, N_k_t_immigration] offset_log_immigration[estimate_immigration, sexes];
  matrix[N_k_a_emigration, N_k_t_emigration] offset_log_emigration[estimate_emigration, sexes];
}
transformed parameters {
  // initialize matrix to store projected population
  matrix<lower = 0>[A, Y + 1] population[sexes];
  real population_input_groups[N_pop];

  // spline offset parameters
  matrix[1, Y] spline_offset_log_srb = rep_matrix(0, 1, Y);
  matrix[A_f, Y] spline_offset_log_asfr = rep_matrix(0, A_f, Y);
  matrix[A, 1] spline_offset_log_baseline[sexes] = rep_array(rep_matrix(0, A, 1), sexes);
  matrix[A + 1, Y] spline_offset_logit_survival[sexes] = rep_array(rep_matrix(0, A + 1, Y), sexes);
  matrix[A, Y] spline_offset_net_migration[sexes] = rep_array(rep_matrix(0, A, Y), sexes);
  matrix[A, Y] spline_offset_log_immigration[sexes] = rep_array(rep_matrix(0, A, Y), sexes);
  matrix[A, Y] spline_offset_log_emigration[sexes] = rep_array(rep_matrix(0, A, Y), sexes);

  // untransformed ccmpp input parameters
  matrix<lower = 0>[1, Y] srb;
  matrix<lower = 0>[A_f, Y] asfr;
  matrix<lower = 0>[A, 1] baseline[sexes];
  matrix<lower = 0, upper = 1>[A + 1, Y] survival[sexes];
  matrix[A, Y] net_migration[sexes];
  matrix<lower = 0>[A, Y] immigration[estimate_immigration, sexes];
  matrix<lower = 0>[A, Y] emigration[estimate_emigration, sexes];

  // calculate spline offsets
  if (estimate_srb) {
    spline_offset_log_srb = offset_log_srb[1] * B_t_srb';
  }
  if (estimate_asfr) {
    spline_offset_log_asfr = B_a_asfr * offset_log_asfr[1] * B_t_asfr';
  }
  for (s in 1:sexes) {
    if (estimate_baseline) {
      spline_offset_log_baseline[s] = B_a_baseline * offset_log_baseline[1, s];
    }
    if (estimate_survival) {
      spline_offset_logit_survival[s] = B_a_survival * offset_logit_survival[1, s] * B_t_survival';
    }
    if (estimate_net_migration) {
      spline_offset_net_migration[s] = B_a_net_migration * offset_net_migration[1, s] * B_t_net_migration';
    }
    if (estimate_immigration) {
      spline_offset_log_immigration[s] = B_a_immigration * offset_log_immigration[1, s] * B_t_immigration';
    }
    if (estimate_emigration) {
      spline_offset_log_emigration[s] = B_a_emigration * offset_log_emigration[1, s] * B_t_emigration';
    }
  }

  // calculate untransformed ccmpp input parameters
  srb = exp(input_log_srb + spline_offset_log_srb);
  asfr = exp(input_log_asfr + spline_offset_log_asfr);
  for (s in 1:sexes) {
    baseline[s] = exp(input_log_baseline[s] + spline_offset_log_baseline[s]);
    survival[s] = inv_logit(input_logit_survival[s] + spline_offset_logit_survival[s]);
    if (estimate_net_migration) {
      net_migration[s] = input_net_migration[s] + spline_offset_net_migration[s];
    } else {
      immigration[estimate_immigration, s] = exp(input_log_immigration[s] + spline_offset_log_immigration[s]);
      emigration[estimate_emigration, s] = exp(input_log_emigration[s] + spline_offset_log_emigration[s]);
      net_migration[s] = immigration[estimate_immigration, s] - emigration[estimate_emigration, s];
    }
  }

  // Level 2 (ccmpp)
  // project population with most detailed age groups
  population = ccmpp(srb, asfr, baseline, survival, net_migration,
                     sexes, interval, A, A_f, A_f_offset, Y);
  // aggregate to census age groups
  population_input_groups = aggregate(population, input_pop_year_index,
                                      input_pop_age_index, input_pop_sex_index,
                                      N_pop);
}
model {
  // LEVEL 4 (informative prior distributions for variance parameters)
  sigma2_population ~ inv_gamma(alpha_population, beta_population);
  if (estimate_srb) {
    sigma2_srb ~ inv_gamma(alpha_srb, beta_srb);
  }
  if (estimate_asfr) {
    sigma2_asfr ~ inv_gamma(alpha_asfr, beta_asfr);
  }
  if (estimate_survival) {
    sigma2_survival ~ inv_gamma(alpha_survival, beta_survival);
  }
  if (estimate_net_migration) {
    sigma2_net_migration ~ inv_gamma(alpha_net_migration, beta_net_migration);
  }
  if (estimate_immigration) {
    sigma2_immigration ~ inv_gamma(alpha_immigration, beta_immigration);
  }
  if (estimate_emigration) {
    sigma2_emigration ~ inv_gamma(alpha_emigration, beta_emigration);
  }
  // LEVEL 3 (model initial estimates ccmpp of ccmpp inputs)

  for (y in 1:N_k_t_srb) {
    if (estimate_srb) {
      offset_log_srb[1, 1, y] ~ normal(0, sqrt(sigma2_srb[1]));
    }
  }
  for (y in 1:N_k_t_asfr) {
    if (estimate_asfr) {
      offset_log_asfr[1, , y] ~ normal(0, sqrt(sigma2_asfr[1]));
    }
  }
  for (s in 1:sexes) {
    if (estimate_baseline) {
      offset_log_baseline[1, s, , 1] ~ normal(0, sqrt(sigma2_population));
    }
  }
  for (y in 1:N_k_t_survival) {
    for (s in 1:sexes) {
      if (estimate_survival) {
        offset_logit_survival[1, s, , y] ~ normal(0, sqrt(sigma2_survival[1]));
      }
    }
  }
  for (y in 1:N_k_t_net_migration) {
    for (s in 1:sexes) {
      if (estimate_net_migration) {
        offset_net_migration[1, s, , y] ~ normal(0, sqrt(sigma2_net_migration[1]));
      }
    }
  }
  for (y in 1:N_k_t_immigration) {
    for (s in 1:sexes) {
      if (estimate_immigration) {
        offset_log_immigration[1, s, , y] ~ normal(0, sqrt(sigma2_immigration[1]));
      }
    }
  }
  for (y in 1:N_k_t_emigration) {
    for (s in 1:sexes) {
      if (estimate_emigration) {
        offset_log_emigration[1, s, , y] ~ normal(0, sqrt(sigma2_emigration[1]));
      }
    }
  }

  // LEVEL 1 (model the census counts)
  input_log_pop_value ~ normal(log(population_input_groups), sqrt(sigma2_population) * input_pop_weight);
}
