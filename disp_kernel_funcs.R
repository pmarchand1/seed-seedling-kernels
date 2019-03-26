# Functions to calculate various quantities related to the dispersal kernels 
#  (PDF, CDF, mode, median and mean dispersal distance)

library(tibble) # for 'lst'

# Probability density by unit area
disp_pdf <- function(kern, r, a, k) {
    switch(kern,
        `2Dt` = k / (pi*a) * (1 + r^2/a)^(-1-k),
        exppow = k / (2*pi*a^2 * gamma(2/k)) * exp(-(r/a)^k),
        invpow = (k-1)*k / (2*pi*a^2) * (1 + r/a)^(-1-k),
        lognorm = 1/(2*pi*r) * dlnorm(r, meanlog = a, sdlog = k),
        weibull = 1/(2*pi*r) * dweibull(r, shape = k, scale = a),
        stop("Undefined kernel")
    )
}

# Cumulative radial probability (0 to r)
disp_cdf <- function(kern, r, a, k) {
    switch(kern,
        `2Dt` = 1 - (1 + r^2/a)^(-k),
        exppow = pgamma((r/a)^k, 2/k),
        invpow = 1 - a^(k-1) * (k*r + a) * (r + a)^(-k),
        lognorm = plnorm(r, meanlog = a, sdlog = k),
        weibull = pweibull(r, shape = k, scale = a),
        stop("Undefined kernel")
    )
}

# Distance from origin to 2D mode
disp_mode <- function(kern, a, k) {
    switch(kern,
        `2Dt` = 0,
        exppow = 0,
        invpow = 0,
        lognorm = exp(a - 2*k^2),
        weibull = ifelse(k > 2, a * ((k-2)/k)^(1/k), 0),
        stop("Undefined kernel")
    )
}

disp_median <- function(kern, a, k) {
    
    # Bounds for numerical approximation
    rmin <- 0.05
    rmax <- 10000
    
    if (kern %in% c("invpow", "exppow")) {
        disp_med <- NA
        cdf_med <- function(r) disp_cdf(kern, r, a, k) - 1/2
        try(disp_med <- unlist(uniroot(cdf_med, c(rmin, rmax))$root), 
            silent = TRUE)
    } else {
        disp_med <- switch(kern,
                       `2Dt` = sqrt(a * (2^(1/k) - 1)),
                       lognorm = exp(a),
                       weibull = a * log(2)^(1/k),           
                       stop("Undefined kernel")
        )
    }
    disp_med
}

disp_mean <- function(kern, a, k) {
    switch(kern,
        `2Dt` = sqrt(pi*a)/2 * gamma(k - 1/2) / gamma(k),
        exppow = a * gamma(3/k) / gamma(2/k),
        invpow = 2*a / (k-2),
        lognorm = exp(a + k^2/2),
        weibull = a * gamma(1 + 1/k),
        stop("Undefined kernel")
    )
}

# This functions returns one draw from the kernel parameters' prior distributions
#  (used to simulate from model)
draw_params <- function(kern, priors_list) {
    if (kern == "2Dt") {
        alpha <- rnorm(1, priors_list$p_alpha[1], priors_list$p_alpha[2])
        inv_k_real <- rnorm(1, priors_list$p_inv_k_real[1], priors_list$p_inv_k_real[2])
        lst(alpha, inv_k_real)
    } else if (kern == "exppow") {
        alpha <- rnorm(1, priors_list$p_alpha[1], priors_list$p_alpha[2])
        inv_k <- rgamma(1, priors_list$p_inv_k[1], priors_list$p_inv_k[2])
        lst(alpha, inv_k)
    } else if (kern == "invpow") {
        alpha <- rnorm(1, priors_list$p_alpha[1], priors_list$p_alpha[2])
        k_real <- rnorm(1, priors_list$p_k_real[1], priors_list$p_k_real[2])
        lst(alpha, k_real)
    } else if (kern == "lognorm") {
        mu_disp <- rnorm(1, priors_list$p_mu_disp[1], priors_list$p_mu_disp[2])
        sd_disp <- rgamma(1, priors_list$p_sd_disp[1], priors_list$p_sd_disp[2])
        lst(mu_disp, sd_disp)
    } else if (kern == "weibull") {
        alpha <- rnorm(1, priors_list$p_alpha[1], priors_list$p_alpha[2])
        inv_k <- rgamma(1, priors_list$p_inv_k[1], priors_list$p_inv_k[2])
        lst(alpha, inv_k)
    } else {
        stop("Unidentified kernel")
    }
}

# This function transforms the parameters on which priors are set (params)
#  to the parameters used in the kernel functions above (a, k)
calc_derived_params <- function(kern, params) {
    if (kern == "2Dt") {
        inv_k <- 2 / (1 + exp(-params$inv_k_real))
        a <- exp(params$alpha)
        k <- 1 / inv_k
        lst(a, k)
    } else if (kern == "exppow") {
        a <- exp(params$alpha - params$inv_k)
        k <- 1 / params$inv_k
        lst(a, k)
    } else if (kern == "invpow") {
        a <- exp(params$alpha)
        k <- 2 + 3 / (1 + exp(-params$k_real))
        lst(a, k)
    } else if (kern == "lognorm") {
        lst(a = params[[1]], k = params[[2]])
    } else if (kern == "weibull") {
        a <- exp(params$alpha - log(log(2)) * params$inv_k)
        k <- 1 / params$inv_k
        lst(a, k)
    } else {
        stop("Unidentified kernel")
    }    
}
