# Calculate various model-averaged quantities based on stacking weights

library(tidyverse)

# Functions used in this script
source("scripts_post_process/disp_curve_funcs.R")

# Load valid models with stacking weights
res_wgt <- read_csv("res_loo_wgts.csv")

# Remove species with lack of valid models at seed or recruit stage
res_wgt <- filter(res_wgt, !str_detect(sp_name, "Lacm|Jaca|Lue"))

# Only keep models that add up to 95% of total weight, nest by species and stage
res_wgt <- group_by(res_wgt, stage, sp_name) %>%
    arrange(sp_name, stage, desc(weight)) %>%
    mutate(prev_cumul_w = coalesce(lag(cumsum(weight)), 0)) %>%
    filter(prev_cumul_w <= 0.95)

res_nest <- select(res_wgt, stage, species = sp_name, disp_mod, err_mod, weight) %>%
    nest(.key = "models")

# Dispersal kernels
avg_disp <- select(res_nest, -models)
avg_disp$disp <- pmap(res_nest, avg_disp_curve,
                      rmin = 1, rmax = 500, npts = 50, logr = TRUE)
avg_disp <- unnest(avg_disp)

saveRDS(avg_disp, "mod_avg_disp.rds")

# Relative survival
sp_list <- set_names(unique(res_nest$species))
surv_df <- map_dfr(sp_list, 
                   ~ avg_rel_surv(filter(res_nest, species == .), 
                                  rmin = 1, rmax = 100, npts = 50, logr = TRUE),
                   .id = "species")
saveRDS(surv_df, "mod_avg_surv.rds")

# Repeat with subset of r values
surv_df2 <- map_dfr(sp_list, 
                    ~ avg_rel_surv(filter(res_nest, species == .), 
                                   rval = c(1, 2, 5, 10, 20, 50, 100)),
                   .id = "species")

saveRDS(surv_df2, "mod_avg_surv_limitedr.rds")


# Dispersal distance summary statistics
avg_stats <- select(res_nest, -models)
avg_stats$stat <- pmap(res_nest, avg_disp_stats)
avg_stats <- unnest(avg_stats) 

saveRDS(avg_stats, "mod_avg_stats.rds")


# Non-dispersal parameters
beta_theta <- mutate(res_nest, 
    nb_w = map_dbl(models, ~ weighted.mean(.$err_mod == "nb", .$weight))) %>%
    select(-models) # nb_w: total weight of neg.binom. count model
beta_theta$params <- pmap(res_nest, avg_beta_theta)
beta_theta <- unnest(beta_theta)

saveRDS(beta_theta, "mod_avg_beta_theta.rds")

