// File: pwm_RP_gene_specific_priors.stan
data {
  int<lower=1> G;                // number of genes
  int<lower=1> N_mu;             // number of mu points (e.g. 9)
  matrix[G, N_mu] R_obs;         // observed RNA abundances (raw, positive)
  matrix[G, N_mu] P_obs;         // observed protein abundances (raw, positive)
  vector[N_mu] mu;               // growth rates
  matrix[G, N_mu] phi_x;         // per gene-phase factor (phi_x[g,i] = phi_hoc(mu_i) or phi_loc)
  vector[N_mu] phi_hoc;          // known phi_hoc(mu_i)

  // gene-specific prior centers and sds (log-space)
  vector[G] mu_logk;
  vector[G] sd_logk;
  //vector[G] max_logk;
  //vector[G] min_logk;

  vector[G] mu_logl;
  vector[G] sd_logl;

  vector[G] mu_logdR;
  vector[G] sd_logdR;

  vector[G] mu_logdP;
  vector[G] sd_logdP;
}
parameters {
  vector[G] logk;    // log(k_g)
  vector[G] logl;    // log(l_g)
  vector[G] logdR;   // log(delta_R_g)
  vector[G] logdP;   // log(delta_P_g)

  real<lower=0> sigma_R; // measurement noise on log scale
  real<lower=0> sigma_P;
}
model {
  // ---- priors (gene-specific) ----
  for (g in 1:G) {
    logk[g]  ~ normal(mu_logk[g], sd_logk[g]); // T[min_logk[g], max_logk[g]];
    logl[g]  ~ normal(mu_logl[g], sd_logl[g]);
    logdR[g] ~ normal(mu_logdR[g], sd_logdR[g]);
    logdP[g] ~ normal(mu_logdP[g], sd_logdP[g]);
  }

  // weak priors for measurement noise
  sigma_R ~ normal(0, 0.5);    // tune as needed
  sigma_P ~ normal(0, 0.5);

  // ---- likelihood (log-space) ----
  for (g in 1:G) {

    real k = exp(logk[g]);
    real l = exp(logl[g]);
    real dR = exp(logdR[g]);
    real dP = exp(logdP[g]);

    for (i in 1:N_mu) {
      // safety floor to avoid log(0)
      real Rpred = phi_x[g, i] * k / (mu[i] + dR);
      real Ppred = Rpred * phi_hoc[i] * l / (mu[i] + dP);

      Rpred = fmax(Rpred, 1e-12);
      Ppred = fmax(Ppred, 1e-12);

      target += normal_lpdf(log(R_obs[g, i] + 1e-12) | log(Rpred), sigma_R);
      target += normal_lpdf(log(P_obs[g, i] + 1e-12) | log(Ppred), sigma_P);
    }
  }
}
generated quantities {
  // optional: return derived means of R and P per gene & mu
  matrix[G, N_mu] R_fit;
  matrix[G, N_mu] P_fit;
  for (g in 1:G) {
    real k = exp(logk[g]);
    real l = exp(logl[g]);
    real dR = exp(logdR[g]);
    real dP = exp(logdP[g]);
    for (i in 1:N_mu) {
      R_fit[g,i] = phi_x[g,i] * k / (mu[i] + dR);
      P_fit[g,i] = R_fit[g,i] * phi_hoc[i] * l / (mu[i] + dP);
    }
  }
}
