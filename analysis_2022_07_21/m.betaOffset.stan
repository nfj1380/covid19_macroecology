// generated with brms 2.17.0
functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  vector<lower=0>[N] denom;  // response denominator
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  // data for spline s(meanAge, k = 6, bs = "ts")
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data for spline s(popDen, k = 6, bs = "ts")
  int nb_2;  // number of bases
  int knots_2[nb_2];  // number of knots
  // basis function matrices
  matrix[N, knots_2[1]] Zs_2_1;
  // data for spline s(rain, k = 6, bs = "ts")
  int nb_3;  // number of bases
  int knots_3[nb_3];  // number of knots
  // basis function matrices
  matrix[N, knots_3[1]] Zs_3_1;
  // data for spline s(temp, k = 6, bs = "ts")
  int nb_4;  // number of bases
  int knots_4[nb_4];  // number of knots
  // basis function matrices
  matrix[N, knots_4[1]] Zs_4_1;
  // data for spline s(PerUrb, k = 6, bs = "ts")
  int nb_5;  // number of bases
  int knots_5[nb_5];  // number of knots
  // basis function matrices
  matrix[N, knots_5[1]] Zs_5_1;
  // data for spline s(HCexpend, k = 6, bs = "ts")
  int nb_6;  // number of bases
  int knots_6[nb_6];  // number of knots
  // basis function matrices
  matrix[N, knots_6[1]] Zs_6_1;
  // data for spline s(Diabetes, k = 6, bs = "ts")
  int nb_7;  // number of bases
  int knots_7[nb_7];  // number of knots
  // basis function matrices
  matrix[N, knots_7[1]] Zs_7_1;
  // data for spline s(cardioDR, k = 6, bs = "ts")
  int nb_8;  // number of bases
  int knots_8[nb_8];  // number of knots
  // basis function matrices
  matrix[N, knots_8[1]] Zs_8_1;
  // data for spline s(spatLag, k = 6, bs = "ts")
  int nb_9;  // number of bases
  int knots_9[nb_9];  // number of knots
  // basis function matrices
  matrix[N, knots_9[1]] Zs_9_1;
  // data for spline s(log_tests, k = 6, bs = "ts")
  int nb_10;  // number of bases
  int knots_10[nb_10];  // number of knots
  // basis function matrices
  matrix[N, knots_10[1]] Zs_10_1;
  // data for spline s(HIV.AIDS, k = 6, bs = "tp")
  int nb_11;  // number of bases
  int knots_11[nb_11];  // number of knots
  // basis function matrices
  matrix[N, knots_11[1]] Zs_11_1;
  // data for spline s(Genital.herpes, k = 6, bs = "tp")
  int nb_12;  // number of bases
  int knots_12[nb_12];  // number of knots
  // basis function matrices
  matrix[N, knots_12[1]] Zs_12_1;
  // data for spline s(Ascariasis, k = 6, bs = "tp")
  int nb_13;  // number of bases
  int knots_13[nb_13];  // number of knots
  // basis function matrices
  matrix[N, knots_13[1]] Zs_13_1;
  // data for spline s(Malaria, k = 6, bs = "tp")
  int nb_14;  // number of bases
  int knots_14[nb_14];  // number of knots
  // basis function matrices
  matrix[N, knots_14[1]] Zs_14_1;
  // data for spline s(Trichuriasis, k = 6, bs = "tp")
  int nb_15;  // number of bases
  int knots_15[nb_15];  // number of knots
  // basis function matrices
  matrix[N, knots_15[1]] Zs_15_1;
  // data for spline s(Tuberculosis, k = 6, bs = "tp")
  int nb_16;  // number of bases
  int knots_16[nb_16];  // number of knots
  // basis function matrices
  matrix[N, knots_16[1]] Zs_16_1;
  // data for spline s(Hookworm, k = 6, bs = "tp")
  int nb_17;  // number of bases
  int knots_17[nb_17];  // number of knots
  // basis function matrices
  matrix[N, knots_17[1]] Zs_17_1;
  // data for spline s(Schistosomiasis, k = 6, bs = "tp")
  int nb_18;  // number of bases
  int knots_18[nb_18];  // number of knots
  // basis function matrices
  matrix[N, knots_18[1]] Zs_18_1;
  // data for spline s(Lymphatic.filariasis, k = 6, bs = "tp")
  int nb_19;  // number of bases
  int knots_19[nb_19];  // number of knots
  // basis function matrices
  matrix[N, knots_19[1]] Zs_19_1;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_shape_2;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  // log response denominator
  vector[N] log_denom = log(denom);
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // spline coefficients
  // parameters for spline s(meanAge, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  // parameters for spline s(popDen, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_2[1]] zs_2_1;
  real<lower=0> sds_2_1;  // standard deviations of spline coefficients
  // parameters for spline s(rain, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_3[1]] zs_3_1;
  real<lower=0> sds_3_1;  // standard deviations of spline coefficients
  // parameters for spline s(temp, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_4[1]] zs_4_1;
  real<lower=0> sds_4_1;  // standard deviations of spline coefficients
  // parameters for spline s(PerUrb, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_5[1]] zs_5_1;
  real<lower=0> sds_5_1;  // standard deviations of spline coefficients
  // parameters for spline s(HCexpend, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_6[1]] zs_6_1;
  real<lower=0> sds_6_1;  // standard deviations of spline coefficients
  // parameters for spline s(Diabetes, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_7[1]] zs_7_1;
  real<lower=0> sds_7_1;  // standard deviations of spline coefficients
  // parameters for spline s(cardioDR, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_8[1]] zs_8_1;
  real<lower=0> sds_8_1;  // standard deviations of spline coefficients
  // parameters for spline s(spatLag, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_9[1]] zs_9_1;
  real<lower=0> sds_9_1;  // standard deviations of spline coefficients
  // parameters for spline s(log_tests, k = 6, bs = "ts")
  // standarized spline coefficients
  vector[knots_10[1]] zs_10_1;
  real<lower=0> sds_10_1;  // standard deviations of spline coefficients
  // parameters for spline s(HIV.AIDS, k = 6, bs = "tp")
  // standarized spline coefficients
  vector[knots_11[1]] zs_11_1;
  real<lower=0> sds_11_1;  // standard deviations of spline coefficients
  // parameters for spline s(Genital.herpes, k = 6, bs = "tp")
  // standarized spline coefficients
  vector[knots_12[1]] zs_12_1;
  real<lower=0> sds_12_1;  // standard deviations of spline coefficients
  // parameters for spline s(Ascariasis, k = 6, bs = "tp")
  // standarized spline coefficients
  vector[knots_13[1]] zs_13_1;
  real<lower=0> sds_13_1;  // standard deviations of spline coefficients
  // parameters for spline s(Malaria, k = 6, bs = "tp")
  // standarized spline coefficients
  vector[knots_14[1]] zs_14_1;
  real<lower=0> sds_14_1;  // standard deviations of spline coefficients
  // parameters for spline s(Trichuriasis, k = 6, bs = "tp")
  // standarized spline coefficients
  vector[knots_15[1]] zs_15_1;
  real<lower=0> sds_15_1;  // standard deviations of spline coefficients
  // parameters for spline s(Tuberculosis, k = 6, bs = "tp")
  // standarized spline coefficients
  vector[knots_16[1]] zs_16_1;
  real<lower=0> sds_16_1;  // standard deviations of spline coefficients
  // parameters for spline s(Hookworm, k = 6, bs = "tp")
  // standarized spline coefficients
  vector[knots_17[1]] zs_17_1;
  real<lower=0> sds_17_1;  // standard deviations of spline coefficients
  // parameters for spline s(Schistosomiasis, k = 6, bs = "tp")
  // standarized spline coefficients
  vector[knots_18[1]] zs_18_1;
  real<lower=0> sds_18_1;  // standard deviations of spline coefficients
  // parameters for spline s(Lymphatic.filariasis, k = 6, bs = "tp")
  // standarized spline coefficients
  vector[knots_19[1]] zs_19_1;
  real<lower=0> sds_19_1;  // standard deviations of spline coefficients
  real Intercept_shape;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  real log_beta ; // slope for offset/rate
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // actual spline coefficients
  vector[knots_2[1]] s_2_1;
  // actual spline coefficients
  vector[knots_3[1]] s_3_1;
  // actual spline coefficients
  vector[knots_4[1]] s_4_1;
  // actual spline coefficients
  vector[knots_5[1]] s_5_1;
  // actual spline coefficients
  vector[knots_6[1]] s_6_1;
  // actual spline coefficients
  vector[knots_7[1]] s_7_1;
  // actual spline coefficients
  vector[knots_8[1]] s_8_1;
  // actual spline coefficients
  vector[knots_9[1]] s_9_1;
  // actual spline coefficients
  vector[knots_10[1]] s_10_1;
  // actual spline coefficients
  vector[knots_11[1]] s_11_1;
  // actual spline coefficients
  vector[knots_12[1]] s_12_1;
  // actual spline coefficients
  vector[knots_13[1]] s_13_1;
  // actual spline coefficients
  vector[knots_14[1]] s_14_1;
  // actual spline coefficients
  vector[knots_15[1]] s_15_1;
  // actual spline coefficients
  vector[knots_16[1]] s_16_1;
  // actual spline coefficients
  vector[knots_17[1]] s_17_1;
  // actual spline coefficients
  vector[knots_18[1]] s_18_1;
  // actual spline coefficients
  vector[knots_19[1]] s_19_1;
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_shape_2;
  real lprior = 0;  // prior contributions to the log posterior
  vector[N] beta_denom;
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
  // compute actual spline coefficients
  s_2_1 = sds_2_1 * zs_2_1;
  // compute actual spline coefficients
  s_3_1 = sds_3_1 * zs_3_1;
  // compute actual spline coefficients
  s_4_1 = sds_4_1 * zs_4_1;
  // compute actual spline coefficients
  s_5_1 = sds_5_1 * zs_5_1;
  // compute actual spline coefficients
  s_6_1 = sds_6_1 * zs_6_1;
  // compute actual spline coefficients
  s_7_1 = sds_7_1 * zs_7_1;
  // compute actual spline coefficients
  s_8_1 = sds_8_1 * zs_8_1;
  // compute actual spline coefficients
  s_9_1 = sds_9_1 * zs_9_1;
  // compute actual spline coefficients
  s_10_1 = sds_10_1 * zs_10_1;
  // compute actual spline coefficients
  s_11_1 = sds_11_1 * zs_11_1;
  // compute actual spline coefficients
  s_12_1 = sds_12_1 * zs_12_1;
  // compute actual spline coefficients
  s_13_1 = sds_13_1 * zs_13_1;
  // compute actual spline coefficients
  s_14_1 = sds_14_1 * zs_14_1;
  // compute actual spline coefficients
  s_15_1 = sds_15_1 * zs_15_1;
  // compute actual spline coefficients
  s_16_1 = sds_16_1 * zs_16_1;
  // compute actual spline coefficients
  s_17_1 = sds_17_1 * zs_17_1;
  // compute actual spline coefficients
  s_18_1 = sds_18_1 * zs_18_1;
  // compute actual spline coefficients
  s_19_1 = sds_19_1 * zs_19_1;
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[, 1];
  r_1_shape_2 = r_1[, 2];
  beta_denom = exp(log_beta)*denom;
  lprior += student_t_lpdf(Intercept | 3, 11.2, 2.6);
  lprior += student_t_lpdf(sds_1_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_2_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_3_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_4_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_5_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_6_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_7_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_8_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_9_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_10_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_11_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_12_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_13_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_14_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_15_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_16_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_17_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_18_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(sds_19_1 | 3, 0, 2.6)
    - 1 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += student_t_lpdf(Intercept_shape | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.6)
    - 2 * student_t_lccdf(0 | 3, 0, 2.6);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
  lprior += normal_lpdf(log_beta|0,5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1 + Zs_4_1 * s_4_1 + Zs_5_1 * s_5_1 + Zs_6_1 * s_6_1 + Zs_7_1 * s_7_1 + Zs_8_1 * s_8_1 + Zs_9_1 * s_9_1 + Zs_10_1 * s_10_1 + Zs_11_1 * s_11_1 + Zs_12_1 * s_12_1 + Zs_13_1 * s_13_1 + Zs_14_1 * s_14_1 + Zs_15_1 * s_15_1 + Zs_16_1 * s_16_1 + Zs_17_1 * s_17_1 + Zs_18_1 * s_18_1 + Zs_19_1 * s_19_1;
    // initialize linear predictor term
    vector[N] shape = Intercept_shape + rep_vector(0.0, N);
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      shape[n] += r_1_shape_2[J_1[n]] * Z_1_shape_2[n];
    }
    for (n in 1:N) {
      // apply the inverse link function
      shape[n] = exp(shape[n]);
    }
    target += neg_binomial_2_log_lpmf(Y | mu + log(beta_denom), shape .* beta_denom);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_2_1);
  target += std_normal_lpdf(zs_3_1);
  target += std_normal_lpdf(zs_4_1);
  target += std_normal_lpdf(zs_5_1);
  target += std_normal_lpdf(zs_6_1);
  target += std_normal_lpdf(zs_7_1);
  target += std_normal_lpdf(zs_8_1);
  target += std_normal_lpdf(zs_9_1);
  target += std_normal_lpdf(zs_10_1);
  target += std_normal_lpdf(zs_11_1);
  target += std_normal_lpdf(zs_12_1);
  target += std_normal_lpdf(zs_13_1);
  target += std_normal_lpdf(zs_14_1);
  target += std_normal_lpdf(zs_15_1);
  target += std_normal_lpdf(zs_16_1);
  target += std_normal_lpdf(zs_17_1);
  target += std_normal_lpdf(zs_18_1);
  target += std_normal_lpdf(zs_19_1);
  target += std_normal_lpdf(to_vector(z_1));
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_shape_Intercept = Intercept_shape;
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
