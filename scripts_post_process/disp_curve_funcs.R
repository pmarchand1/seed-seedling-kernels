# Functions to calculate model-averaged dispersal kernels, relative survival,
#  dispersal distance statistics, and non-dispersal parameters
#  from a list of candidate models with weights

library(tidyverse)

source("disp_kernel_funcs.R")

# - For a given lifestage, species, kernel (disp_mod) and count_model(err_mod),
#   calculate the dispersal kernels for a vector of distance values (rval),
#   OR create this vector from rmin, rmax, npts.
#   (if logr = TRUE, use geometric series from rmin to rmax)
# - For each posterior sample and r value, calculates three probability types:
#    2D prob. density (dens2D),
#    2D prob. density times mean production per basal area (densmu),
#    radial cumulative probability (cumul)
# - Outputs a data frame with one row by sample, r value and prob. type
make_disp_curve <- function(stage, species, disp_mod, err_mod, rval = NULL,
                            rmin = 1, rmax = 50, npts = 50, logr = FALSE) {
    # Load result file
    res <- readRDS(paste0("results/", stage, "_", species, "_", disp_mod,
                          "_", err_mod, ".rds"))

    nsamp <- length(res$samples[[1]])
    
    disp_pars <- res$samples[1:2] # First two params in samples are for dispersal model
    deriv_pars <- calc_derived_params(disp_mod, disp_pars)
    pars_df <- cbind(data.frame(samp = 1:nsamp), deriv_pars,
                     res$samples[c("mu_beta", "sd_beta")]) %>%
        map_df(as.vector) # Remove attributes from columns
    
    # Values of r to calculate dispersal at
    # (if rval vector not supplied, create from rmin, rmax, npts)
    if (is.null(rval)) {
        if (logr) {
            rval <- exp(seq(log(rmin), log(rmax), length.out = npts))
            
        } else {
            rval <- seq(rmin, rmax, length.out = npts)    
        }    
    } else {
        npts <- length(rval)
    }
    
    # Create data frame to calcualte probability for each combination of r and pars
    disp_df <- pars_df[rep(1:nsamp, each = npts), ]
    disp_df$r <- rep(rval, nrow(pars_df))
    
    disp_df <- mutate(disp_df,
                      dens2D = disp_pdf(disp_mod, r, a, k),
                      densmu = dens2D * exp(mu_beta + sd_beta^2/2),
                      cumul = disp_cdf(disp_mod, r, a, k))

    disp_df <- gather(disp_df, key = "prob_type", value = "prob", 
                      dens2D, densmu, cumul)
}

# - For a given lifestage, species, and a data frame "models" 
#   with columns disp_mod, err_mod and weight, calculates a model-averaged
#   version of the output of make_disp_curve by resampling from the samples of
#   each model proportionally to weights.
# - Arguments rval, rmin, rmax, npts and logr are passed to make_disp_curve.
# - If return_summary = TRUE, returns a data frame with quantiles
#   (2.5%, 25%, 50%, 75%, 97.5%) for each r value and prob. type
avg_disp_curve <- function(stage, species, models, rval = NULL,
                           rmin = 1, rmax = 50, npts = 50, logr = FALSE,
                           return_summary = TRUE) {
    # Get dispersal curve samples for all candidate models
    nmod <- nrow(models)
    models <- mutate(models, 
        disp_df = map2(disp_mod, err_mod, ~ make_disp_curve(stage, species, .x, .y, 
                                                            rval, rmin, rmax, npts, logr))) %>%
        unnest()
    
    # Resample from combined sample proportionally to stacking weights
    models <- select(models, prob_type, r, weight, prob) %>%
        group_by(prob_type, r) %>%
        nest() %>%
        mutate(resamp_df = map(data, ~ data.frame(
            isamp = seq_len(nrow(.) / nmod),
            resamp = sample(.$prob, size = nrow(.) / nmod, replace = TRUE, prob = .$weight)
        ))) %>%
        select(-data) %>% 
        unnest()
    
    if (return_summary) {
        # Return quantiles of averaged posterior predictive distribution
        models <- group_by(models, prob_type, r) %>%
            summarize(lo95 = quantile(resamp, 0.025), lo50 = quantile(resamp, 0.25), 
                      med = quantile(resamp, 0.5), hi50 = quantile(resamp, 0.75), 
                      hi95 = quantile(resamp, 0.975))        
    }
    models
}

# - Calculates the model-averaged relative survival from seed to recruit
#   and (if applicable) recruit to seedling for a given species.
# - Input res_df is a data frame with columns stage, species (constant),
#   and models (list-column where each row contains a data frame),
#   which are passed as arguments to avg_disp_curve, along with the
#   additional arguments rval to logr.
# - After getting the combined sample for each lifestage, this function calculates
#   the ratio of recruit/seed and seedling/recruit densities for each sample,
#   then outputs summary quantiles for each stage transition and r value.
# - Relative survival values are normalized to mean of 0 on a log scale
#   (i.e. over the values of r considered, median relative survival is set to 1)
avg_rel_surv <- function(res_df, rval = NULL,
                         rmin = 1, rmax = 50, npts = 50, logr = FALSE) {
    # Get dispersal curve samples for all candidate models
    avg_disp <- pmap(res_df, avg_disp_curve, rval, rmin, rmax, npts, logr,
                     return_summary = FALSE)
    avg_disp <- map(avg_disp, ~ filter(., prob_type == "dens2D")) %>%
        set_names(res_df$stage) %>%
        bind_rows(.id = "stage") %>%
        select(-prob_type)
    # Calculate relative survival based on 2D density ratios
    rel_surv <- spread(avg_disp, key = stage, value = resamp) %>%
        mutate(`seed-to-recruit` = log(recruit) - log(seed))
    if ("seedling" %in% res_df$stage) {
        rel_surv <- mutate(rel_surv, 
                           `recruit-to-seedling` = log(seedling) - log(recruit)) %>%
            select(-seedling)
    }
    rel_surv <- select(rel_surv, -isamp, -seed, -recruit) %>%
        gather(-r, key = "stage", value = "log_surv")
    # Normalize log(rel.surv.) to mean of 0, then calculate quantiles
    rel_surv <- group_by(rel_surv, stage) %>%
        mutate(log_rel_surv = log_surv - mean(log_surv)) %>%
        mutate(rel_surv = exp(log_rel_surv)) %>%
        group_by(stage, r) %>%
        summarize(lo95 = quantile(rel_surv, 0.025), lo50 = quantile(rel_surv, 0.25), 
                  med = quantile(rel_surv, 0.5), hi50 = quantile(rel_surv, 0.75), 
                  hi95 = quantile(rel_surv, 0.975)) %>%
        ungroup()
}

# For a given lifestage, species and model (disp_mod, err_mod),
#  calculates modal, median and mean dispersal distance for each posterior sample
# Note: modal distance is the distance corresponding to the 2D density mode
get_disp_stats <- function(stage, species, disp_mod, err_mod) {

    res <- readRDS(paste0("results/", stage, "_", species, "_", disp_mod,
                          "_", err_mod, ".rds"))
    
    if ("rmed" %in% names(res$samples)) {
        disp_stats <- data.frame(res$samples[c("rmode", "rmean", "rmed")])
    } else {
        disp_pars <- res$samples[1:2]
        deriv_pars <- calc_derived_params(disp_mod, disp_pars)
        rmed <- pmap_dbl(deriv_pars, ~ disp_median(disp_mod, ..1, ..2))
        disp_stats <- data.frame(res$samples[c("rmode", "rmean")], rmed)
    }
    disp_stats <- map_df(disp_stats, as.vector) # Remove attributes from columns
    
    disp_stats <- gather(disp_stats, key = "stat", value = "value", rmode, rmed, rmean)
}

# Averages the output of get_disp_stats for multiple models
#  and returns summary quantiles for each statistic
avg_disp_stats <- function(stage, species, models) {

    # Get dispersal summary statistics for all candidate models
    models <- models %>%
        mutate(disp_df = map2(disp_mod, err_mod, 
                              ~ get_disp_stats(stage, species, .x, .y))) %>%
        unnest()
    
    # Resample from combined sample proportionally to stacking weights
    models <- select(models, stat, weight, value) %>%
        group_by(stat) %>%
        nest() %>%
        mutate(resamp = map(data, ~ sample(.$value, size = nrow(.), 
                                           replace = TRUE, prob = .$weight))) %>%
        select(-data) %>% 
        unnest()
    
    # Return quantiles of averaged posterior predictive distribution
    group_by(models, stat) %>%
        summarize(lo95 = quantile(resamp, 0.025), lo50 = quantile(resamp, 0.25), 
                  med = quantile(resamp, 0.5), hi50 = quantile(resamp, 0.75), 
                  hi95 = quantile(resamp, 0.975))
}

# For a given lifestage, species and model (disp_mod, err_mod),
#  extracts posterior samples of the non-dispersal parameters
get_beta_theta <- function(stage, species, disp_mod, err_mod) {
    res <- readRDS(paste0("results/", stage, "_", species, "_", disp_mod,
                          "_", err_mod, ".rds"))
    res <- res$samples[names(res$samples) %in% c("mu_beta", "sd_beta", "ri_theta")] %>%
        map_df(as.vector)
    # If clumping parameter absent (Poisson), it is equivalent to 0
    if (!("ri_theta" %in% colnames(res))) 
        res$ri_theta <- 0

    gather(res, key = "param", value = "value")
}

# Obtain model-averaged summary quantiles for the non-dispersal parameters
avg_beta_theta <- function(stage, species, models) {

    models <- models %>%
        mutate(pars = map2(disp_mod, err_mod, 
                           ~ get_beta_theta(stage, species, .x, .y))) %>%
        unnest()
    
    # Resample from combined sample proportionally to stacking weights
    models <- select(models, param, weight, value) %>%
        group_by(param) %>%
        nest() %>%
        mutate(resamp = map(data, ~ sample(.$value, size = nrow(.), 
                                           replace = TRUE, prob = .$weight))) %>%
        select(-data) %>% 
        unnest()
    
    # Return quantiles of averaged posterior predictive distribution
    group_by(models, param) %>%
        summarize(lo95 = quantile(resamp, 0.025), lo50 = quantile(resamp, 0.25), 
                  med = quantile(resamp, 0.5), hi50 = quantile(resamp, 0.75), 
                  hi95 = quantile(resamp, 0.975))
}
