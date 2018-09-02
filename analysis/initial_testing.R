library(tidyverse)
library(magrittr)
library(parallel)
library(readxl)
library(rstan)

cdat = read_excel('data/nm.4219-S2.xlsx')

# HPAF, ASPC1, PATU are cell lines
# TXX is time point XX
# other differences are replicates

cdat %<>%
  select(matches('GENE'), HPAF_T0, HPAF_T35A, HPAF_T35B) # just look at one time point for now

fit_nb = function(counts){
  counts = counts$count
  fn_to_min = function(param_vec){
    # param_vec[1] nb mean
    # param_vec[2] nb size
    -sum(dnbinom(counts,
                 mu = param_vec[1],
                 size = param_vec[2],
                 log = TRUE))
  }

  stats::nlminb(start = c(100, 1),
                objective = fn_to_min,
                lower = rep(.Machine$double.xmin, 2))
}

fit_gamma = function(param_estimates){

  fn_to_min = function(ab_vec){
    -sum(dgamma(param_estimates,
                shape = ab_vec[1],
                rate = ab_vec[2],
                log = TRUE))
  }

  stats::nlminb(start = c(1,1),
                objective = fn_to_min,
                lower = rep(.Machine$double.xmin, 2))
}

sample_depths = cdat %>%
  select(matches('HPAF')) %>%
  map_df(~sum(.x) / 1e6) %>%
  gather(sample, depth_factor)

nb_fits = cdat %>%
  gather(sample, count, matches('HPAF')) %>%
  group_by(GENE, sample) %>%
  nest(.key = counts) %>%
  mutate(n_gRNA = map_int(counts, nrow)) %>%
  filter(n_gRNA >= 3) %>%
  mutate(nb_fit = mclapply(counts, fit_nb, mc.cores = 3)) %>%
  mutate(converged = map_lgl(nb_fit, ~.x$convergence == 0),
         mu_est = map_dbl(nb_fit, ~.x$par[1]),
         phi_est = map_dbl(nb_fit, ~.x$par[2])) %>%
  left_join(sample_depths, by = 'sample') %>%
  mutate(depth_adj_mu_est = mu_est / depth_factor)


#### We'll use T0 oas the prior
gamma_priors = nb_fits %>%
  filter(grepl('T0|T35A|T35B', sample)) %>%
  mutate(differing_values = map_lgl(counts, ~n_distinct(.x$count) > 2)) %>% # This cuts out 17 / 15657
  filter(differing_values) %>%
  select(-nb_fit, -counts) %>%
  filter(converged) %>%
  group_by(sample) %>%
  nest(.key = nb_params) %>%
  mutate(mu_gamma_prior = mclapply(nb_params, function(.x){fit_gamma(.x$depth_adj_mu_est)}, mc.cores = 3),
         phi_gamma_prior = mclapply(nb_params, function(.x){fit_gamma(.x$phi_est)}, mc.cores = 3),
         mu_converged = map_lgl(mu_gamma_prior, ~.x$convergence == 0),
         phi_converged = map_lgl(phi_gamma_prior, ~.x$convergence == 0))


make_prior_line = function(prior_list,
                           nb_params,
                           param = 'mu'){

  param_range = nb_params %>%
    select(depth_adj_mu_est) %>%
    .[,1] %>%
    range

  if (param == 'mu') {
    data_frame(x = seq(param_range[1],
                       param_range[2],
                       length.out = 1000),
               y = dgamma(x,
                          shape = prior_list$par[1],
                          rate = prior_list$par[2]))
  } else {
    data_frame(x = 10**seq(log10(param_range[1]),
                       log10(param_range[2]),
                       length.out = 1000),
               y = dgamma(x,
                          shape = prior_list$par[1],
                          rate = prior_list$par[2]))
  }

}
prior_lines = gamma_priors %>%
  mutate(mu_prior_line = map2(mu_gamma_prior, nb_params, make_prior_line),
         phi_prior_line = map2(phi_gamma_prior, nb_params, make_prior_line, param = 'phi')) %>%
  select(sample, matches('line')) %>%
  select_all(~gsub('_prior_line', '', .)) %>%
  gather(param, prior_dens, -sample) %>%
  unnest

prior_lines %<>%
  bind_rows(prior_lines %>% mutate(sample = rep('HPAF_T35A', n())),
            prior_lines %>% mutate(sample = rep('HPAF_T35B', n())))

nb_fits %>%
  filter(converged) %>%
  ggplot(aes(depth_adj_mu_est)) +
  geom_histogram(aes(y = ..density..),
                 bins = 50) +
  facet_wrap('sample',
             scales = 'free_x') +
  geom_line(data = prior_lines %>% filter(param == 'mu') %>% rename(mu_est = x),
            aes(mu_est, y)) +
  labs(title = 'CRISPR Screen MLE NB Mean parameters with fitted Gamma priors',
       subtitle = 't0 prior reused for other time points')

nb_fits %>%
  filter(converged) %>%
  mutate(differing_values = map_lgl(counts, ~n_distinct(.x$count) > 2)) %>%
  filter(differing_values) %>%
  ggplot(aes(phi_est)) +
  geom_histogram(aes(y = ..density..),
                 bins = 50) +
  facet_wrap('sample',
             scales = 'free') +
  geom_line(data = prior_lines %>% filter(param == 'phi') %>% rename(phi_est = x),
            aes(phi_est, y)) +
  labs(title = 'CRISPR Screen MLE NB Phi parameters with Gamma priors') +
  scale_x_log10()

#### Stan model ----
# needs to have per gRNA effects? not sure how input concentration is accounted for

crispr_model = stan_model('src/stan_files/crispr_model.stan')

analyze_gene = function(gene_dat){
  data_list = list(n_gRNA = nrow(gene_dat),
                   treatment_counts = gene_dat %>% pull(HPAF_T35B) %>% as.integer,
                   ctrl_counts = gene_dat %>% pull(HPAF_T35A) %>% as.integer,
                   t0_counts = gene_dat %>% pull(HPAF_T0) %>% as.integer,
                   ctrl_depth = sample_depths$depth_factor[2],
                   treatment_depth = sample_depths$depth_factor[3],
                   t0_depth = sample_depths$depth_factor[1],
                   t0_mu_prior = gamma_priors %>% filter(sample == 'HPAF_T0') %>% .$mu_gamma_prior %>% .[[1]] %>% .$par,
                   t0_phi_prior =  gamma_priors %>% filter(sample == 'HPAF_T0') %>% .$phi_gamma_prior %>% .[[1]] %>% .$par,
                   treatment_mu_prior =  gamma_priors %>% filter(sample == 'HPAF_T35B') %>% .$mu_gamma_prior %>% .[[1]] %>% .$par,
                   treatment_phi_prior =  gamma_priors %>% filter(sample == 'HPAF_T35B') %>% .$phi_gamma_prior %>% .[[1]] %>% .$par,
                   ctrl_mu_prior =  gamma_priors %>% filter(sample == 'HPAF_T35A') %>% .$mu_gamma_prior %>% .[[1]] %>% .$par,
                   ctrl_phi_prior =  gamma_priors %>% filter(sample == 'HPAF_T35A') %>% .$phi_gamma_prior %>% .[[1]] %>% .$par)

  samp_res = sampling(object = crispr_model,
           data = data_list,
           chains = 1,
           iter = 10000,
           warmup = 500)

  coda::HPDinterval(coda::mcmc(rstan::extract(samp_res, 'log_fc')$log_fc %>% as.matrix(ncol = 1)))
}

lfc_hdis = cdat %>%
  group_by(GENE) %>%
  nest %>%
  head(n = 30) %>%
  mutate(log_fc_hdi = mclapply(data, analyze_gene, mc.cores = 3))
