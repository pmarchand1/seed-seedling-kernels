# This script runs one instance of the hierarchical Bayesian model
#  (suitable for parallelization on a computing cluster)
#
# Call: Rscript run_single_model.R lifestage species disp_mod err_mod
# - lifestage: seed, recruit or seedling
# - species: scientific name (w/ underscore betweeen words)
# - disp_mod: dispersal kernel (2Dt, exppow, invpow, lognorm or weibull)
# - err_mod: count model (nb or pois)

library(methods)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(rstan)
library(bayesplot) # for 'rhat' function

source("data_subset_funcs.R") # Load data and subsetting functions
source("calc_dists_weights_func.R") # Functions for edge-correction weights

# Get command-line arguments
argv <- commandArgs(trailingOnly=TRUE)
lifestage <- argv[1]
species <- gsub("_", " ", argv[2])
disp_mod <- argv[3]
err_mod <- argv[4]
# lifestage <- "recruit"
# species <- "Beilschmiedia tovarensis"
# disp_mod <- "weibull"
# err_mod <- "nb"
print(paste(lifestage, species, disp_mod, err_mod))

# Fixed parameters
dist_trap_max <- 10 # For seedlings, max. distance from a seed trap
trap_area <- switch(lifestage, seed = 0.5, recruit = 1, seedling = 1)
years <- switch(lifestage, seed = 1988:2014, recruit = 1994:2014, 
                seedling = seedling_years)

# Load priors and arrange into list
disp_priors <- read_csv("disp_priors.csv") %>%
    filter(model == disp_mod) %>%
    select(-model)
repr_priors <- read_csv("repr_priors.csv") %>%
    filter(stage == lifestage) %>%
    select(-stage)
priors <- bind_rows(disp_priors, repr_priors)
if (err_mod != "nb") 
    priors <- filter(priors, param != "p_ri_theta")
priors_list <- setNames(map2(priors$prior_mean, priors$prior_sd, c), 
                        priors$param)

# Iterations and number of chains for Stan
n_warmup <- 600
n_iter <- 1600 # includes warmup
n_chains <- 2

# Extract tree data, create tree size x year matrix
tree_data <- map_dfr(set_names(years), ~ subset_trees(species, .), .id = "year")
tree_data <- spread(tree_data, key = "year", value = "size", fill = 0)
tree_size <- as.matrix(select(tree_data, -id, -x, -y))

# Extract seed, recruit or seedling data
subset_func <- switch(lifestage,
    seed = subset_seeds,
    recruit = subset_recruits,
    seedling = subset_seedlings
)
seed_data <- map_dfr(set_names(years), ~ subset_func(species, .), .id = "year")

# Reshape seed_data differently for seedlings vs. seeds or recruits
if (lifestage == "seedling") {
    seed_data <- spread(seed_data, key = "year", value = "seeds")
    # Only keep within dist_trap_max of a seed trap
    seed_traps <- read_csv("data/seed_traps.csv")
    dist_trap <- sqrt(outer(seed_data$x, seed_traps$X, "-")^2 + 
                          outer(seed_data$y, seed_traps$Y, "-")^2)
    seed_data <- seed_data[apply(dist_trap, 1, min) < dist_trap_max, ]
    seed_data$itrap <- 1:nrow(seed_data)
    # Scale by dt to get seedlings by year
    nseed <- as.matrix(select(seed_data, -itrap, -x, -y))
    nseed <- scale(nseed, center = FALSE, scale = seedling_dt)
} else { 
    # Seed and recruit data is kept in "sparse matrix" format 
    #   due to change in number of traps over years
    seed_data <- mutate(seed_data, iyear = as.integer(as.factor(year)), 
                                   itrap = as.integer(as.factor(trap)))
    nseed <- seed_data$seeds    
}

# Round fractions up or down randomly (with p = fraction)
round_frac <- function(x) {
    ifelse(runif(length(x)) < (x %% 1), ceiling(x), floor(x))
}
nseed[nseed %% 1 > 0] <- round_frac(nseed[nseed %% 1 > 0])

# Calculate tree-trap distance matrix, maximum radius in plot and edge-correction weights
traps <- distinct(seed_data, itrap, x, y)
dist_weights <- calc_dists_weights(traps, tree_data)

# Other input variables for model
total_size <- colSums(tree_size)
plot_area <- (xmax - xmin) * (ymax - ymin)
size_density <- total_size / plot_area

# Create data list for Stan
if (lifestage == "seedling") {
    data_list <- lst(trap_area, nyear = ncol(tree_size), ntree = nrow(tree_size), 
                     ntrap = nrow(traps), tree_size, size_density, 
                     nseed = as.vector(nseed))        
} else {
    data_list <- lst(trap_area, nyear = ncol(tree_size), ntree = nrow(tree_size), 
                     ntrap = nrow(traps), tree_size, size_density, 
                     iyear = seed_data$iyear, itrap = seed_data$itrap, 
                     nobs = length(nseed), nseed)    
}
data_list <- c(data_list, dist_weights, priors_list)

# Check for missing data
if (any(is.na(unlist(data_list)))) stop("Missing values in data.")

# Run Stan model
model_file <- paste("stan_models/disp", disp_mod, err_mod, 
                    ifelse(lifestage == "seedling", "mat.stan", "sparse.stan"), 
                    sep = "_")
res <- stan(model_file, data = data_list, chains = n_chains, 
            warmup = n_warmup, iter = n_iter, cores = 2)

# Export diagnostics, LOO results and parameter samples
pars_keep <- c("alpha|inv_k|k_real|mu_disp|sd_disp|mu_beta|sd_beta|ri_theta")
diags <- list(ndiv = get_num_divergent(res),
              max_tree = get_num_max_treedepth(res),
              bfmi = min(get_bfmi(res)),
              max_rhat = max(rhat(res, regex_pars = pars_keep)))

loo_res <- loo(res)

samples <- extract(res)
samples <- samples[!(names(samples) == "log_lik")]

res_sum <- list(samples = samples, diags = diags, loo = loo_res)

saveRDS(res_sum, paste0("results/", lifestage, "_", species, "_", 
                        disp_mod, "_", err_mod, ".rds"))