# Generate prior distribution for median and mean distances by kernel

library(rstan)

mod_list <- c("2Dt", "exppow", "invpow", "lognorm", "weibull")

model_files <- paste0("stan_models/test_prior_", mod_list,".stan")

# Produce data list for each model
disp_priors <- read_csv("disp_priors.csv") %>%
    nest(-model)
priors_list <- map(disp_priors$data, 
                   ~ setNames(map2(.$prior_mean, .$prior_sd, c), .$param))

# Iterate over model files and data lists for each of 5 models, run in Stan
prpr <- map2(model_files, priors_list,
             ~ stan(file = .x, data = .y, chains = 1, warmup = 0, 
                    iter = 10000, refresh = 10000, algorithm = "Fixed_param"))

# Extract samples and combine in one data frame
pars <- map(prpr, extract)
names(pars) <- mod_list

rstats <- map_dfr(pars, ~ .[c("rmed", "rmean")], .id = "model")

rstats <- gather(rstats, -model, key = "stat", value = "value") %>%
    mutate(model = recode(model, `2Dt` = "2D t", exppow = "Exp. power", 
                          invpow = "Inv. power", lognorm = "Log-normal",
                          weibull = "Weibull"),
           stat = recode(stat, rmed = "median", rmean = "mean"))
rstats$stat <- factor(rstats$stat, levels = c("median", "mean"))

# Plot prior distributions
ggplot(rstats, aes(x = value)) +
    labs(x = "r (m)", y = "Prior probability density", color = "", fill = "") +
    geom_density(aes(fill = stat, color = stat), alpha = 0.5) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_log10() +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.3, 3000)) +
    facet_wrap(~ model) +
    theme(axis.title.y = element_text(margin = margin(r = 20)),
          legend.position = c(0.75, 0.25),
          strip.background = element_blank(), 
          strip.text = element_text(margin = margin(b = 10)),
          panel.spacing = unit(14, "pt"))

group_by(rstats, stat, model) %>%
    summarize(med = median(value), 
              q01 = quantile(value, 0.01, na.rm = TRUE),
              q99 = quantile(value, 0.99, na.rm = TRUE))

# 1%: 0.6 to 1.5 for med, 0.8 to 2.4 for mean
# 99%: 70 to 300 for med, 200 to 600 for mean


