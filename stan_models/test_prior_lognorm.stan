data {
    real p_mu_disp[2];
    real p_sd_disp[2];
}

generated quantities {
    
    real mu_disp;
    real<lower=0> sd_disp;
    real rmed;
    real rmean;
    
    mu_disp = normal_rng(p_mu_disp[1], p_mu_disp[2]);
    sd_disp = gamma_rng(p_sd_disp[1], p_sd_disp[2]);
    
    rmed = exp(mu_disp);
    rmean = exp(mu_disp + sd_disp^2/2);
}
