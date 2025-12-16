//   Copyright 2022 James Watson and Phrutsamon Wongnak, Mahidol Oxford Tropical Medicine Research Unit
//
//   Licensed under a CC-BY License
/**
Model 0 used for the simulation studies with no covariate adjustment and no RNaseP
**/



data {
  // Patient data
  int<lower=0> Ntot;                           // Number of PCR data points
  int<lower=0,upper=Ntot> N_obs;               // Number of PCR data points with viral load less than LOD
  int<lower=0> n_id;                           // Number of individuals
  int<lower=1,upper=n_id> id[Ntot];            // Patient identifier for each PCR sample
  int<lower=1,upper=Ntot> ind_start[n_id];     // Starting index for each patient: used in generated quantities block
  real<lower=0> Time_max;                      // Max follow-up time
  real<lower=0,upper=Time_max> obs_day[Ntot];  // Time since randomisation for sample
  real log_10_vl[Ntot];                        // log base 10 viral load in copies per mL
  real log10_cens_vl[Ntot];                    // censoring value for censored observation
  int<lower=1> K_trt;                          // Number of treatment arms
  matrix[Ntot,K_trt] trt_mat;                  // Trt matrix

  // priors
  real alpha_0_prior_mean; // prior mean intercept
  real alpha_0_prior_sd;   // prior sd intercept

  real beta_0_prior_mean;  // prior mean slope
  real beta_0_prior_sd;    // prior sd slope

  real trt_effect_sd;   // prior sd on treatment effect

  real sigma_logvl_mean;   // prior mean of sd in error model
  real sigma_logvl_sd;     // prior sd of sd in error model

}

transformed data {
  vector[2] zeros2;
  for(i in 1:2) zeros2[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[2] L_Omega;     // correlation matrix
  vector<lower=0>[2] sigmasq_u;        // variance of random effects

  // Measurement error
  real<lower=0> sigma_logvl;

  // Population parameters
  real alpha_0;                             // population intercept
  real beta_0;                              // population slope
  vector[K_trt-1] trt_effect;               // Estimates of the treatment effect

  // Random effects
  vector[2] theta_rand_id[n_id];            // individual random effects vector

  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}

transformed parameters {
  real pred_log10_vl[Ntot];
  vector[Ntot] trt_slope;
  {
    vector[K_trt] trt_effect_prime;

    trt_effect_prime = append_row(0, trt_effect);
    trt_slope = trt_mat * trt_effect_prime;

    // calculate predicted log viral load under the model parameters
    for(i in 1:Ntot){
      pred_log10_vl[i] =
      alpha_0 + theta_rand_id[id[i]][1] +
      beta_0*exp(trt_slope[i]+theta_rand_id[id[i]][2])*obs_day[i];
    }
  }
}

model {
  //***** Prior *****
  // error model - degrees of freedom
  t_dof ~ exponential(1);
  // error model variance
  sigma_logvl ~ normal(sigma_logvl_mean,sigma_logvl_sd) T[0,];

  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2); // covariance matrix - random effects for individs
  // individual random effects
  for(i in 1:n_id) theta_rand_id[i] ~ multi_normal_cholesky(zeros2, diag_pre_multiply(sigmasq_u, L_Omega));

  // Population parameters
  alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  beta_0 ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  trt_effect ~ normal(0,trt_effect_sd); // Treatment effect

  //***** Likelihood *****
  // Non censored observations
  log_10_vl[1:N_obs] ~ student_t(t_dof, pred_log10_vl[1:N_obs], sigma_logvl);

  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }

}
