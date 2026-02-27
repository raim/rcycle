data {
  int<lower=1> G;               // number of genes
  int<lower=1> N_mu;            // number of mu values (e.g., 9)
  matrix[G, N_mu] R_obs;        // observed RNA abundances
  matrix[G, N_mu] P_obs;        // observed protein abundances
  vector[N_mu] mu;              // growth rates
  matrix[G, N_mu] phi_x;        // phase factor for each gene and mu
  vector[N_mu] phi_hoc;         // HOC phase factor

}

parameters {
  // Hierarchical parameters (log-space)

  real mu_logk;
  real<lower=0> tau_logk;
  vector[G] z_logk;

  real mu_logdR;
  real<lower=0> tau_logdR;
  vector[G] z_logdR;

  real mu_logl;
  real<lower=0> tau_logl;
  vector[G] z_logl;

  real mu_logdP;
  real<lower=0> tau_logdP;
  vector[G] z_logdP;

  // empirical prior distributions
  real logk_mu;
  real logk_sd;
  real logdr_mu;
  real logdr_sd;
  real logl_mu;
  real logl_sd;
  real logdp_mu;
  real logdp_sd;

  real<lower=0> sigma_R;       // RNA measurement error
  real<lower=0> sigma_P;       // protein measurement error
}

transformed parameters {
  vector[G] logk = mu_logk + tau_logk * z_logk;
  vector[G] logdR = mu_logdR + tau_logdR * z_logdR;
  vector[G] logl = mu_logl + tau_logl * z_logl;
  vector[G] logdP = mu_logdP + tau_logdP * z_logdP;
}

model {

  // Priors
  mu_logk ~ normal(logk_mu, logk_sd);
  mu_logdR ~ normal(logdr_mu, logdr_sd);
  mu_logl ~ normal(logl_mu, logl_sd);
  mu_logdP ~ normal(logdp_mu, logdp_sd);

  tau_logk ~ cauchy(0, 1);
  tau_logdR ~ cauchy(0, 1);
  z_logk ~ normal(0, 1);
  z_logdR ~ normal(0, 1);
  tau_logl ~ cauchy(0, 1);
  tau_logdP ~ cauchy(0, 1);
  z_logl ~ normal(0, 1);
  z_logdP ~ normal(0, 1);
  
  sigma_R ~ normal(0, 0.5);
  sigma_P ~ normal(0, 0.5);

  // Likelihood
  for (g in 1:G) {

    real k = exp(logk[g]);
    real dR = exp(logdR[g]);
    real l = exp(logl[g]);
    real dP = exp(logdP[g]);

    for (i in 1:N_mu) {
      real R_pred = phi_x[g, i] * k / (mu[i] + dR);
      real P_pred = R_pred * phi_hoc[i] * l / (mu[i] + dP);

      target += normal_lpdf(log(R_obs[g, i]+1e-8) | log(R_pred), sigma_R);
      target += normal_lpdf(log(P_obs[g, i]+1e-8) | log(P_pred), sigma_P);
    }
  }
}
