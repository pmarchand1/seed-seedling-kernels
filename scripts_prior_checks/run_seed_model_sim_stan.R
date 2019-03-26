# This script runs one simulation of seed dispersal model from priors,
#  and fits model to simulted data with Stan
#  (suitable for parallelization on a computing cluster)
#
# Call: Rscript run_seed_model_sim_stan.R species disp_mod err_mod sim_id
# - species: 6-letter abbreviation of species name from BCI list
# - disp_mod: dispersal kernel (2Dt, exppow, invpow, lognorm or weibull)
# - err_mod: count model (nb or pois)
# - sim_id: simulation ID (to save results from multiple sims from same model)

library(methods)
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(rstan)
library(bayesplot) # for 'rhat' function


# Prepare inputs for simulation -------------------------------------------

# Load data and required functions
source("load_sim_data.R")
source("calc_dists_weights_funcs.R")
source("disp_kernel_funcs.R")

# Fixed parameters
trap_area <- 0.5
nyear <- 5

# Get command-line arguments
argv <- commandArgs(trailingOnly=TRUE)
species <- argv[1]
disp_mod <- argv[2]
err_mod <- argv[3]
sim_id <- argv[4]
# species <- "BEILPE"
# disp_mod <- "weibull"
# err_mod <- "nb"
print(paste(species, disp_mod, err_mod))

# Load priors and arrange into list
disp_priors <- read_csv("disp_priors.csv") %>%
    filter(model == disp_mod) %>%
    select(-model)
repr_priors <- read_csv("repr_priors.csv")
priors <- bind_rows(disp_priors, repr_priors)
if (err_mod != "nb") 
    priors <- filter(priors, param != "p_ri_theta")
priors_list <- setNames(map2(priors$prior_mean, priors$prior_sd, c), 
                        priors$param)

# Iterations and number of chains for Stan
n_warmup <- 500
n_iter <- 1000 # includes warmup
n_chains <- 2

# Get trees from species
trees <- filter(trees, sp_code6 == species)
tree_size <- trees$rba_cm2

# Calculate tree-trap distance matrix, maximum radius in plot and edge-correction weights
dist_weights <- calc_dists_weights(traps, trees)

# Other input variables for model
total_size <- sum(tree_size)
plot_area <- (xmax - xmin) * (ymax - ymin)
size_density <- total_size / plot_area


# Simulate seed counts from model -----------------------------------------

ntree <- nrow(trees)
ntrap <- nrow(traps)

mu_beta <- rnorm(1, priors_list$p_mu_beta[1], priors_list$p_mu_beta[2])
sd_beta <- abs(rnorm(1, priors_list$p_sd_beta[1], priors_list$p_sd_beta[2]))

params <- draw_params(disp_mod, priors_list)
derived_pars <- calc_derived_params(disp_mod, params)

disp_kern <- dist_weights$wgt * disp_pdf(disp_mod, dist_weights$r, 
                                         derived_pars[[1]], derived_pars[[2]])

offplot <- 1 - disp_cdf(disp_mod, dist_weights$rmax, 
                        derived_pars[[1]], derived_pars[[2]])

beta_off <- rnorm(nyear, 0, 1)
b <- exp(mu_beta + sd_beta * beta_off)

mu <- trap_area * (as.vector(disp_kern %*% tree_size) + offplot * size_density) %o% b

if (err_mod == "nb") {
    ri_theta <- abs(rnorm(1, priors_list$p_ri_theta[1], priors_list$p_ri_theta[2]))
    theta <- ri_theta^(-2)
    nseed <- rnbinom(length(mu), size = theta, mu = as.vector(mu))
} else {
    nseed <- rpois(length(mu), lambda = as.vector(mu))
}

# Infer parameters in Stan ------------------------------------------------

data_list <- lst(trap_area, nyear, ntree, ntrap, tree_size, size_density, nseed)
data_list <- c(data_list, dist_weights, priors_list)

# Check for missing data
if (any(is.na(unlist(data_list)))) stop("Missing values in data.")

model_file <- paste("stan_models/disp", disp_mod, err_mod, "fixed.stan", 
                    sep = "_")
res <- stan(model_file, data = data_list, chains = n_chains, 
            warmup = n_warmup, iter = n_iter, cores = 1)

# Save original parmeters, simulated seed counts, 
#  posterior samples and HMC diagnostics

if (err_mod == "nb") {
    pars <- unlist(c(params, lst(mu_beta, sd_beta, ri_theta)))
} else {
    pars <- unlist(c(params, lst(mu_beta, sd_beta)))
}

rhat_vals <- possibly(function() {
    rh <- rhat(res)
    rh <- rh[names(pars)]
}, NA)()

diags <- list(ndiv = get_num_divergent(res),
              max_tree = get_num_max_treedepth(res),
              bfmi = min(get_bfmi(res)),
              rhat = max(rhat_vals))

samples <- possibly(function() {
    samp <- extract(res)
    samp <- samp[!(names(samp) == "log_lik")]
}, NULL)()

saveRDS(lst(pars, nseed, diags, samples), 
        paste0("sim_res/", species, "_", disp_mod, "_", err_mod, sim_id, ".rds"))
