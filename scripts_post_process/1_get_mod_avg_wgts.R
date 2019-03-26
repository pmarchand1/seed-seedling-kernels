# From Stan model outputs, calculate stacking weights for each species and lifestage

library(tidyverse)
library(loo)
library(coda)

# Function to calculate stacking weights from a data frame of model results
#  with "loo" column corresponding to output of loo function
calc_stack_wgts <- function(mod_df) {
    if (nrow(mod_df) == 1) {
        mod_df$weight <- 1
    } else {
        mod_df$weight <- invoke(cbind, map(mod_df$loo, ~ .$pointwise[, 1])) %>%
            stacking_weights()        
    }
    select(mod_df, -loo)
}

# List result files and extract model info.
res_df <- tibble(file = dir("results")) %>%
    separate(file, into = c("stage", "sp_name", "disp_mod", "err_mod", "ext"),
             sep = "[_\\.]", remove = FALSE) %>%
    select(-ext)

# Extract total counts and # of non-zero counts by species and lifestage
#  to compare with posterior predictive values
sp_tab <- read_csv("sp_tab.csv")
sp_counts <- select(sp_tab, sp_name, seedling_ntot = seedling_near_traps, 
                    seedling_nnz = seedling_near_trap_yrs,
                    recruit_ntot = recruit_count, recruit_nnz = recruit_trap_yrs,
                    seed_ntot = seed_count, seed_nnz = seed_trap_yrs) %>%
    gather(-sp_name, key = "stage_var", value = "count") %>%
    separate(stage_var, into = c("stage", "var"), sep = "_") %>%
    spread(key = var, value = count)

res_df <- inner_join(res_df, sp_counts)

# Load all result files in list
all_res <- map(file.path("results", res_df$file), readRDS)

# Add Stan diagnostics, goodness-of-fit tests and loo output
res_df <- bind_cols(res_df, map_dfr(all_res, "diags"))
res_df <- mutate(res_df, 
    ll_pval = map_dbl(all_res, ~ mean(.$samples$pval)),
    nnz_pval = map2_dbl(all_res, nnz, ~ mean(.x$samples$nnz < .y)),
    ntot_pval = map2_dbl(all_res, ntot, ~ mean(.x$samples$tot < .y))
)
res_df$loo <- map(all_res, "loo")

# Define valid models as those with no more than one divergence or max tree depth,
#  bfmi > 0.2, and non-extreme goodness-of-fit p-values
res_df <- mutate(res_df, 
    valid = ndiv <= 1 & max_tree <= 1 & bfmi >= 0.2 & 
            ll_pval >= 0.05 & nnz_pval >= 0.025 & nnz_pval <= 0.975 & 
            ntot_pval >= 0.025 & ntot_pval <= 0.975)

# Select valid models and nest data frame by species and lifestage
res_valid <- filter(res_df, valid) %>%
    nest(-stage, -sp_name)

res_wgt <- suppressWarnings(
    mutate(res_valid, data = map(data, calc_stack_wgts)) %>%
    unnest()
)

write_csv(res_wgt, "res_loo_wgts.csv")

# Model comparison table for appendix
mod_tab <- left_join(res_df, select(res_wgt, file, weight)) %>%
    select(sp_name, stage, disp_kernel = disp_mod, count_mod = err_mod, 
           n_div = ndiv, n_maxtree = max_tree, bfmi, ll_pval, nnz_pval, 
           ntot_pval, max_rhat, valid, weight) %>%
    mutate(stage = factor(stage, levels = c("seed", "recruit", "seedling"))) %>%
    arrange(sp_name, stage, -weight)

write_csv(mod_tab, "appendix_model_table.csv")
