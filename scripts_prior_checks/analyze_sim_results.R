# Check consistency and precision of parameter estimates from simulation results

library(tidyverse)
library(rstan)

source("disp_kernel_funcs.R")


# Combine all results -----------------------------------------------------

# This function extracts summary statistics from a set of simulations of the same model
summarize_sims <- function(species, disp_mod, err_mod) {
    # Read sim result files
    res_files <- list.files("sim_res", full.names = TRUE,
                            pattern = paste(species, disp_mod, err_mod, sep = "_"))
    sim_res <- map(res_files, readRDS)
    
    # Extract input parameters and calculate mean dispersal distance
    pars <- map_df(sim_res, ~ as.list(.$pars))
    derived_pars <- calc_derived_params(disp_mod, pars[, 1:2])
    pars$rmean <- map2_dbl(derived_pars[[1]], derived_pars[[2]],
                           ~ disp_mean(disp_mod, .x, .y))
    pars$lrmean <- log(pars$rmean)
    
    # Extract parameter estimates and summarize (posterior mean and s.d.)
    par_est <- map_dfr(sim_res, ~ data.frame(.$samples), .id = "sim") %>%
        mutate(sim = as.integer(sim)) %>%
        mutate(lrmean = log(rmean))
    
    par_sum <- select(par_est, sim, one_of(names(pars))) %>%
        gather(-sim, key = "param", value = "value") %>%
        group_by(param, sim) %>%
        summarize(post_mean = mean(value, na.rm = TRUE), 
                  post_sd = sd(value, na.rm = TRUE))
    
    # Calculate z-values and shrinkage by parameter and simulation
    pars_tall <- mutate(pars, sim = 1:nrow(pars)) %>%
        gather(-sim, key = "param", value = "true_val") %>%
        group_by(param) %>%
        mutate(prior_sd = sd(true_val))
    
    par_sum <- inner_join(pars_tall, par_sum) %>%
        mutate(z_val = (post_mean - true_val) / post_sd,
               shrink = 1 - (post_sd / prior_sd)^2)
    
    # Get diagnostic values and join with parameter estimates
    diag_df <- map_dfr(sim_res, ~ data.frame(.$diags), .id = "sim") %>%
        mutate(sim = as.integer(sim))
    
    res_df <- inner_join(par_sum, diag_df)
}


sim_df <- tibble(species = c(rep("BEILPE", 7), rep("BROSAL", 3)),
                 disp_mod = c("2Dt", "exppow", "invpow", "lognorm", "weibull",
                              "2Dt", "lognorm", "2Dt", "weibull", "2Dt"),
                 err_mod = c(rep("nb", 5), rep("pois", 2), "nb", "nb", "pois"))

sim_df$res <- pmap(sim_df, summarize_sims)
sim_df <- unnest(sim_df)

diag_df <- distinct(sim_df, species, disp_mod, err_mod, sim, ndiv, max_tree, bfmi, rhat)

# Model failure by species and model
fail_df <- group_by(diag_df, species, disp_mod, err_mod) %>%
    summarize(sum(ndiv > 1), sum(max_tree > 1), sum(ndiv > 1 | max_tree > 1))



# Produce graphs ----------------------------------------------------------

# Recode parameter and kernel names for plotting
sim_df$par_rename <- recode(sim_df$param, alpha = "Disp. scale", mu_disp = "Disp. scale",
                            inv_k = "Disp. shape", inv_k_real = "Disp. shape", 
                            k_real = "Disp. shape", sd_disp = "Disp. shape",
                            ri_theta = "Clumping", rmean = "Disp. mean", 
                            lrmean = "Disp. mean", mu_beta = "Fecund. mean",
                            sd_beta = "Fecund. std. dev.")
sim_df$kernel <- recode(sim_df$disp_mod, `2Dt` = "2D t", exppow = "Exp. power",
                        invpow = "Inv. power", lognorm = "Log-normal", weibull = "Weibull")


# Q-Q plot of z-values
ggplot(filter(sim_df, ndiv <= 1, max_tree <= 1, species == "BEILPE", err_mod == "nb",
              param != "lrmean"), aes(sample = z_val, color = kernel)) +
    geom_qq(size = 1, alpha = 0.7) +
    geom_abline(linetype = "dotted") +
    scale_color_brewer(palette = "Dark2") +
    facet_wrap(~ par_rename) +
    theme_bw() +
    theme(strip.background = element_blank())


# Distribution of shrinkage
ggplot(filter(sim_df, ndiv <= 1, max_tree <= 1, param != "rmean", 
              species == "BEILPE", err_mod == "nb"), 
       aes(x = shrink, color = kernel)) +
    labs(x = "Shrinkage", y = "Density") +
    geom_density(bw = 0.1) +
    facet_wrap(~ par_rename, scales = "free_y") +
    coord_cartesian(xlim = c(0, 1)) +
    theme_bw() +
    theme(strip.background = element_blank())

ggplot(filter(sim_df, ndiv <= 1, max_tree <= 1), aes(x = shrink)) +
    stat_ecdf(geom = "point", size = 1) +
    facet_wrap(~ param) +
    scale_x_continuous(limits = c(0.5, 1)) +
    theme_bw()


# Shrinkage vs. true values

ggplot(filter(sim_df, ndiv <= 1, max_tree <= 1, param != "rmean",
              species == "BEILPE", err_mod == "nb"), 
       aes(x = true_val, y = shrink, color = kernel)) +
    labs(x = "Parameter value", y = "Shrinkage") +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~ par_rename, scales = "free_x") +
    theme_bw() +
    theme(strip.background = element_blank())

ggplot(filter(sim_df, ndiv <= 1, max_tree <= 1, param == "lrmean"), 
       aes(x = true_val, y = shrink)) +
    geom_point(size = 1, alpha = 0.3) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~ disp_mod, scales = "free_x") +
    theme_bw()


# Divergences vs. true values

sim_pars <- filter(sim_df, param != "rmean", species == "BEILPE") %>%
    select(species, kernel, err_mod, sim, ndiv, max_tree, par_rename, true_val)
sim_pars <- spread(sim_pars, key = par_rename, value = true_val)

# Effect of kernel parameters for neg.binom. counts
ggplot(filter(sim_pars, err_mod == "nb"), 
       aes(x = `Disp. scale`, y = `Disp. shape`, color = ndiv > 1)) +
    labs(x = "Kernel scale", y = "Kernel shape", color = "> 1 divergent") +
    geom_point(size = 1) +
    facet_wrap(~ kernel, scales = "free") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = c(0.8, 0.2))

# Effect of reproduction parameters for Poisson counts
ggplot(filter(sim_pars, err_mod == "pois"), 
       aes(x = `Fecund. mean`, y = `Fecund. std. dev.`, color = max_tree > 1)) +
    labs(x = expression(mu[beta]), y = expression(sigma[beta]), 
         color = "> 1 max. tree depth") +    
    geom_point(size = 1) +
    facet_wrap(~ kernel, scales = "free") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    theme(strip.background = element_blank())

