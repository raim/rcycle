


data {
  int<lower=1> G;               // number of genes
  int<lower=1> N_mu;            // number of mu values (e.g., 9)
  matrix[G, N_mu] R_obs;        // observed RNA abundances
  matrix[G, N_mu] P_obs;        // observed protein abundances
  vector[N_mu] mu;              // growth rates
  vector[N_mu] phi_hoc;         // HOC phase factor
}

parameters {

  // empirical prior distributions
  real logl_mu;
  real logl_sd;
  real logl_max;
  real logdp_mu;
  real logdp_sd;
  real logdp_max;

  // Hierarchical parameters (log-space)
  real<lower=0, upper=logl_max>  mu_logl;
  real<lower=0> tau_logl;
  vector[G] z_logl;

  real<lower=0, upper=logdp_max> mu_logdP;
  real<lower=0> tau_logdP;
  vector[G] z_logdP;

  real<lower=0> sigma_P;       // protein measurement error
}

transformed parameters {
  vector[G] logl = mu_logl + tau_logl * z_logl;
  vector[G] logdP = mu_logdP + tau_logdP * z_logdP;
}

model {

  // Priors
  mu_logl ~ normal(logl_mu, logl_sd);
  mu_logdP ~ normal(logdp_mu, logdp_sd);

  tau_logl ~ cauchy(0, 1);
  tau_logdP ~ cauchy(0, 1);
  z_logl ~ normal(0, 1);
  z_logdP ~ normal(0, 1);

  sigma_P ~ normal(0, 0.5);

  // Likelihood
  for (g in 1:G) {

    real l = exp(logl[g]);
    real dP = exp(logdP[g]);

    for (i in 1:N_mu) {
      real P_pred = R_obs[g, i] * phi_hoc[i] * l / (mu[i] + dP);

      target += normal_lpdf(log(P_obs[g, i]+1e-8) | log(P_pred), sigma_P);
    }
  }
}
