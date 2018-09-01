data{
  int<lower=0> n_gRNA;
  int<lower=0> treatment_counts[n_gRNA];
  int<lower=0> ctrl_counts[n_gRNA];
  int<lower=0> t0_counts[n_gRNA];
  int<lower=0> ctrl_depth;
  int<lower=0> treatment_depth;
  int<lower=0> t0_depth;

  real<lower=0> t0_mu_prior[2];
  real<lower=0> t0_phi_prior[2];

  // here I allow for the post-screen mu priors to be different,
  // but in practice I think I'm going to use the same priors as the t0 mu prior
  real<lower=0> treatment_mu_prior[2];
  real<lower=0> ctrl_mu_prior[2];

  real<lower=0> treatment_phi_prior[2];
  real<lower=0> ctrl_phi_prior[2];
}
parameters{
  real<lower=0> t0_mu[n_gRNA];
  real<lower=0> t0_phi;

  real<lower=0> ctrl_mu;
  real<lower=0> treatment_mu;
  real<lower=0> ctrl_phi;
  real<lower=0> treatment_phi;
}
model{

  for (g in 1:n_gRNA) {
    t0_mu[g] ~ gamma(t0_mu_prior[1], t0_mu_prior[2]);
  }

  t0_phi ~ gamma(t0_phi_prior[1], t0_phi_prior[2]);

  treatment_mu ~ gamma(treatment_mu_prior[1], treatment_mu_prior[2]);
  treatment_phi ~ gamma(treatment_phi_prior[1], treatment_phi_prior[2]);
  ctrl_mu ~ gamma(ctrl_mu_prior[1], ctrl_mu_prior[2]);
  ctrl_phi ~ gamma(ctrl_phi_prior[1], ctrl_phi_prior[2]);

  for (g in 1:n_gRNA) {
    t0_counts[g] ~ neg_binomial_2(t0_mu[g] * t0_depth, t0_phi);

    treatment_counts[g] ~ neg_binomial_2(t0_mu[g] * treatment_depth * treatment_mu, treatment_phi);
    ctrl_counts[g] ~ neg_binomial_2(t0_mu[g] * ctrl_depth * ctrl_mu, ctrl_phi);
  }
}
generated quantities{
  real log_fc;
  log_fc = log(treatment_mu) - log(ctrl_mu);
}
