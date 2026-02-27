data {
  int<lower=1> G;               // number of genes
  int<lower=1> N_mu;            // number of mu values (e.g., 9)
  matrix[G, N_mu] R_obs;        // observed RNA abundances
  vector[N_mu] mu;              // growth rates
  matrix[G, N_mu] phi_x;        // transcription duty-cycle for each gene
}

parameters {

  // empirical prior distributions
  real logk_mu;
  real logk_sd;
  real logk_max;
  real logk_min;
  real logdr_mu;
  real logdr_sd;
  real logdr_max;
  real logdr_min;

  // Hierarchical parameters (log-space)
  real<lower=logk_min, upper=logk_max> mu_logk;
  real<lower=0> tau_logk;
  vector[G] z_logk;

  //real mu_logdR;
  real<lower=logdr_min, upper=logdr_max> mu_logdR;
  real<lower=0> tau_logdR;
  vector[G] z_logdR;


  real<lower=0> sigma_R;       // RNA measurement error
}

transformed parameters {
  vector[G] logk = mu_logk + tau_logk * z_logk;
  vector[G] logdR = mu_logdR + tau_logdR * z_logdR;
}

model {

  // Priors
  mu_logk ~ normal(logk_mu, logk_sd);
  mu_logdR ~ normal(logdr_mu, logdr_sd);

  tau_logk ~ cauchy(0, 1);
  tau_logdR ~ cauchy(0, 1);
  z_logk ~ normal(0, 1);
  z_logdR ~ normal(0, 1);
  sigma_R ~ normal(0, 0.5);

  // Likelihood
  for (g in 1:G) {
    real k = exp(logk[g]);
    real dR = exp(logdR[g]);

    for (i in 1:N_mu) {
      real R_pred = phi_x[g, i] * k / (mu[i] + dR);

      target += normal_lpdf(log(R_obs[g, i]+1e-8) | log(R_pred), sigma_R);
    }
  }
}
