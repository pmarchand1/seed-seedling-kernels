# Creates all graphs in paper from model-averaged results

library(tidyverse)
library(cowplot)

# Settings and helper functions -------------------------------------------

stage_order <- c("seed", "recruit", "seedling")
stage_trans_order <- c("seed-to-recruit", "recruit-to-seedling")

# Default color scales
scale_colour_discrete <- function(...) scale_colour_brewer(palette = "Dark2", ...)
scale_fill_discrete <- function(...) scale_fill_brewer(palette = "Dark2", ...)

# Abbreviate species name by replacing genus with initial
abbr_genus <- function(sp_name) {
    str_replace(sp_name, "(^[A-Z])[a-z]+", "\\1\\.")
}

# Load additional species data --------------------------------------------

# Seed mass
seed_dat <- read_csv("bci_raw_data/seeds/SeedMassforTraits_20150225.csv")
seed_dat <- select(seed_dat, sp_code4 = SP4, seed_mass = SEED_DRY)
# Strength of distance-dependence from Murphy et al. 2017
murphy <- read_csv("bci_raw_data/murphy2017_tabS1.csv")
murphy <- select(murphy, sp_code6, ndd_0_50)

sp_tab <- read_csv("sp_tab.csv")
sp_tab <- inner_join(sp_tab, seed_dat) %>%
    left_join(murphy) %>%
    select(sp_code6, species = sp_name, tree_count, seed_mass, ndd_0_50)


# Dispersal kernels -------------------------------------------------------

avg_disp <- readRDS("mod_avg_disp.rds")
avg_disp <- filter(avg_disp, prob_type == "dens2D") %>%
    mutate(stage = factor(stage, levels = stage_order),
           species = abbr_genus(species))

pdf("figures/post_kernel.pdf", width = 6, height = 7)
ggplot(avg_disp) +
    geom_line(aes(x = r, y = med, color = stage)) +
    geom_ribbon(aes(x = r, ymin = lo95, ymax = hi95, fill = stage), alpha = 0.3) +
    labs(x = "r (m)", y = expression("Probability density" * (1/m^2)), 
         color = "", fill = "") +
    scale_y_log10(breaks = 10^(seq(-5, -1, 2))) +
    coord_cartesian(ylim = c(max(min(avg_disp$lo95), 1E-6), max(avg_disp$hi95)),
                    xlim = c(1, 200)) +
    facet_wrap(~ species, ncol = 4) +
    theme_bw() +
    theme(strip.background = element_blank(), legend.position = "top",
          axis.text = element_text(color = "black"),
          panel.border = element_rect(color = "grey"),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)))
dev.off()


# Relative survival -------------------------------------------------------

surv_df <- readRDS("mod_avg_surv.rds")

surv_df <- mutate(surv_df, 
    stage = factor(stage, levels = stage_trans_order), 
    species = abbr_genus(species))

pdf("figures/relative_survival.pdf", width = 6, height = 7)
ggplot(surv_df, aes(x = r)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_line(aes(y = med, color = stage)) +
    geom_ribbon(aes(x = r, ymin = lo50, ymax = hi50, fill = stage), alpha = 0.3) +
    labs(x = "r (m)", y = "Relative survival", color = "", fill = "") +
    scale_x_log10(breaks = c(1, 3, 10, 30, 100)) +
    scale_y_log10() +
    facet_wrap(~ species, ncol = 4, scales = "free_y") +
    theme_bw() +
    theme(strip.background = element_blank(), legend.position = "top",          
          axis.text = element_text(color = "black"),
          panel.border = element_rect(color = "grey"),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)))
dev.off()


# Relative survival (discrete) --------------------------------------------
# (compare with negative distance-dependence estimates from Murphy et al. 2017)
surv_discr <- readRDS("mod_avg_surv_limitedr.rds")

surv_discr <- select(surv_discr, species, stage, r, med) %>%
    spread(key = r, value = med) %>%
    mutate(ratio5_50 = `50` /`5`) %>%
    select(species, stage, ratio5_50)

surv_murphy <- inner_join(surv_discr, sp_tab) %>%
    filter(!is.na(ndd_0_50)) %>%
    mutate(stage = factor(stage, levels = stage_trans_order), 
           species = abbr_genus(species))

ggplot(surv_murphy, aes(x = ratio5_50, y = ndd_0_50, color = stage, fill = stage)) +
    labs(x = "Relative survival at 50 m vs. 5 m",
         y = "Seedling negative distance-dependence", color = "", fill = "") +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_x_log10()


# Median and mean distance ------------------------------------------------

avg_stats <- readRDS("mod_avg_stats.rds")

avg_stats <- mutate(avg_stats,
    stat = factor(recode(stat, rmode = "Mode", rmean = "Mean", rmed = "Median"),
                  levels = c("Mode", "Median", "Mean")),
    stage = factor(stage, levels = stage_order),
    species = fct_rev(abbr_genus(species))
)

pdf("figures/mean_med_disp.pdf", width = 6, height = 4.5)
ggplot(filter(avg_stats, stat != "Mode"), aes(x = species, y = med, color = stage)) +
    labs(x = "", y = "Mean or median dispersal distance (m)") +
    geom_point(size = 1.5, position = position_dodge(width = 0.3)) +
    geom_linerange(aes(ymin = lo50, ymax = hi50), size = 1, alpha = 0.5,
                   position = position_dodge(width = 0.3)) +
    geom_linerange(aes(ymin = lo95, ymax = hi95), size = 0, alpha = 0.5,
                   position = position_dodge(width = 0.3)) +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
    facet_wrap(~stat) +
    coord_flip(ylim = c(0.5, 12000)) +
    theme_cowplot(font_size = 10) +
    theme(panel.grid.major.x = element_line(colour = "grey90"),
          strip.background = element_blank(), 
          strip.text = element_text(margin = margin(b = 10)))
dev.off()


# Beta and theta ----------------------------------------------------------

beta_theta <- readRDS("mod_avg_beta_theta.rds")

beta_theta <- mutate(beta_theta,
    stage = factor(stage, levels = stage_order),
    species = fct_rev(species)
)

ggplot(filter(beta_theta, param == "ri_theta"), 
       aes(x = species, y = med, ymin = lo95, ymax = hi95, color = stage)) +
    labs(y = expression("Clumping (" * sqrt(1/theta) * ")"), x = "", color = "") +
    geom_point(size = 2, position = position_dodge(width = 0.3)) +
    geom_linerange(aes(ymin = lo50, ymax = hi50), size = 2, alpha = 0.5,
                   position = position_dodge(width = 0.3)) +
    geom_linerange(aes(ymin = lo95, ymax = hi95), alpha = 0.5,
                   position = position_dodge(width = 0.3)) +
    coord_flip()

ggplot(filter(beta_theta, param == "mu_beta"), 
       aes(x = species, y = med, ymin = lo95, ymax = hi95, color = stage)) +
    labs(y = expression(mu[beta]), x = "", color = "") +
    geom_point(size = 2, position = position_dodge(width = 0.3)) +
    geom_linerange(aes(ymin = lo50, ymax = hi50), size = 2, alpha = 0.5,
                   position = position_dodge(width = 0.3)) +
    geom_linerange(aes(ymin = lo95, ymax = hi95), alpha = 0.5,
                   position = position_dodge(width = 0.3)) +
    coord_flip()

ggplot(filter(beta_theta, param == "sd_beta"), 
       aes(x = species, y = med, ymin = lo95, ymax = hi95, color = stage)) +
    labs(y = expression(sigma[beta]), x = "", color = "") +
    geom_point(size = 2, position = position_dodge(width = 0.3)) +
    geom_linerange(aes(ymin = lo50, ymax = hi50), size = 2, alpha = 0.5,
                   position = position_dodge(width = 0.3)) +
    geom_linerange(aes(ymin = lo95, ymax = hi95), alpha = 0.5,
                   position = position_dodge(width = 0.3)) +
    coord_flip()


# Survival, relative survival and clumping vs. tree count and seed mass ----

# Overall survival rate estimated as ratio of mean production per basal area
#  between consecutive life stages
surv <- filter(beta_theta, param != "ri_theta") %>%
    select(species, stage, param, med) %>%
    spread(key = param, value = med) %>%
    # beta is the log of mean production (log-normal distribution)
    mutate(beta = mu_beta + 0.5 * sd_beta^2) %>%
    select(species, stage, beta) %>%
    spread(key = stage, value = beta) %>%
    mutate(`seed-to-recruit` = recruit - seed,
           `recruit-to-seedling` = seedling - recruit) %>%
    select(-seed, -recruit, -seedling) %>%
    gather(key = "stage", value = "log_surv", -species) %>%
    mutate(surv_prob = exp(log_surv))

# Ratio of clumping factor between consecutive stages
theta_comp <- filter(beta_theta, param == "ri_theta") %>%
    select(species, stage, med) %>%
    spread(key = stage, value = med) %>%
    mutate(`seed-to-recruit` = recruit / seed,
           `recruit-to-seedling` = seedling / recruit) %>%
    select(-seed, -recruit, -seedling) %>%
    gather(key = "stage", value = "clump_ratio", -species)

# Combine with relative survival at 50 m vs. 5 m, and add covariates
surv_cov <- inner_join(surv, surv_discr) %>%
    inner_join(theta_comp) %>%
    inner_join(sp_tab) %>%
    select(species, stage, surv_prob, ratio5_50, clump_ratio, tree_count, seed_mass) %>%
    mutate(stage = factor(stage, levels = stage_trans_order))

# Clumping ratio vs. relative survival
ggplot(surv_cov, aes(x = clump_ratio, y = ratio5_50, color = stage, fill = stage)) +
    labs(y = "Relative survival at 50 m vs. 5 m", x = "Clumping ratio", color = "", fill = "") +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_y_log10()


# Reshape data for matrix plot (3 response variables, 2 covariates)

surv_cov <- surv_cov %>%
    gather(surv_prob, ratio5_50, clump_ratio, key = "resp_var", value = "resp_val") %>%
    gather(tree_count, seed_mass, key = "cov_var", value = "cov_val") %>%
    mutate(cov_var = recode(cov_var, tree_count = "Reproductive trees", 
                            seed_mass = "Seed mass (g)"),
           resp_var = recode(resp_var, surv_prob = "Survival rate",
                             ratio5_50 = "Survival at 50 m vs. 5 m",
                             clump_ratio = "Clumping ratio"))

# Identify significant correlations
surv_cov <- mutate(surv_cov, 
                   signif = stage == "seed-to-recruit" & cov_var == "Seed mass (g)")

scale_num <- scales::format_format(scientific = FALSE, drop0trailing = TRUE)

p1 <- ggplot(filter(surv_cov, resp_var != "Clumping ratio"), 
             aes(x = cov_val, y = resp_val, color = stage, fill = stage)) +
    labs(x = "", y = "", color = "", fill = "") +
    geom_point(size = 1) +
    geom_smooth(method = "lm", aes(linetype = signif), size = 0.5) +
    scale_x_log10(labels = scale_num) +
    scale_y_log10(labels = scale_num) +
    scale_linetype_manual(values = c("dotted", "solid"), guide = "none") +
    facet_grid(resp_var ~ cov_var, scales = "free", as.table = FALSE, switch = "both") +
    theme_bw() +
    theme(legend.position = "top", strip.background = element_blank(),
          strip.placement = "outside", strip.text.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(b = 0, r = 10),
          axis.text = element_text(color = "black"),
          panel.border = element_rect(color = "grey"))

p2 <- ggplot(filter(surv_cov, resp_var == "Clumping ratio"), 
             aes(x = cov_val, y = resp_val, color = stage, fill = stage)) +
    labs(x = "", y = "", color = "", fill = "") +
    geom_point(size = 1) +
    geom_smooth(method = "lm", aes(linetype = signif), size = 0.5) +
    scale_x_log10(labels = scale_num) +
    scale_linetype_manual(values = c("dotted", "solid"), guide = "none") +
    coord_cartesian(ylim = c(0, 3)) +
    facet_grid(resp_var ~ cov_var, scales = "free", as.table = FALSE, switch = "both") +
    theme_bw() +
    theme(legend.position = "none", strip.background = element_blank(),
          strip.placement = "outside", plot.margin = margin(t = 0, r = 10),
          axis.text = element_text(color = "black"),
          panel.border = element_rect(color = "grey"))

pdf("figures/trees_seedmass.pdf", width = 6, height = 4.5)
plot_grid(p1, p2, ncol = 1, rel_heights = c(2, 1), align = "v")
dev.off()
