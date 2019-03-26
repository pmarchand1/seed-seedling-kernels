// Seed dispersal model with lognormal kernel

functions {
    //  Calculate expected number of seeds per trap and year
    matrix calc_mu(real mu_disp, real sd_disp, row_vector beta_off, real mu_beta, 
                   real sd_beta, data real trap_area, int nyear, int ntree, int ntrap, 
                   data matrix r, data matrix log_wgt_2pr, data vector rmax, 
                   data matrix tree_size, data row_vector size_density) {
        row_vector[nyear] b; // Fecundity parameter by year
        matrix[ntrap, ntree] disp_kern; // Dispersal kernel
        vector[ntrap] offplot; // Radial integral of disp. kern from rmax to inf.
        matrix[ntrap, nyear] mu;
        
        for (j in 1:ntree) {
            for (i in 1:ntrap) {
                disp_kern[i, j] = exp(log_wgt_2pr[i, j] + 
                                      lognormal_lpdf(r[i, j] | mu_disp, sd_disp));
            }
        }
    
        for (i in 1:ntrap) {
            offplot[i] = 1 - lognormal_cdf(rmax[i], mu_disp, sd_disp);    
        }
    
        b = exp(mu_beta + sd_beta * beta_off);
        
        mu = trap_area * rep_matrix(b, ntrap) .* (disp_kern * tree_size + 
                                                        offplot * size_density);
        return(mu);
    }
}

data {
    // Area of seed traps (meters)
    real<lower=0> trap_area; 
    
    // Number of observations, years, trees and traps
    int<lower=1> nyear;
    int<lower=1> ntree;
    int<lower=1> ntrap;
    
    // Trap-tree dist. matrix, and weights for offplot area
    matrix<lower=0>[ntrap, ntree] r;
    matrix<lower=0>[ntrap, ntree] wgt;
    
    // Maximum distance from each trap to edge of plot
    vector<lower=0>[ntrap] rmax;
    
    // Size of trees, and total tree size density (sum of sizes / plot area)
    matrix<lower=0>[ntree, nyear] tree_size;
    row_vector<lower=0>[nyear] size_density;
    
    // Number of seeds in traps
    int<lower=0> nseed[ntrap*nyear];
    
    // Hyperparameters for parameter priors
    real p_mu_disp[2];
    real p_sd_disp[2];
    real p_mu_beta[2];
    real p_sd_beta[2];
}

transformed data {
    matrix[ntrap, ntree] log_wgt_2pr = log(wgt) - (log2() + log(pi()) + log(r));
}

parameters {
    real mu_disp; // Mean log of kernel
    real sd_disp; // Std.dev. of log of kernel
    row_vector[nyear] beta_off; // Fecundity offset by year
    real mu_beta; // Mean log of b
    real<lower=0> sd_beta; // Std. dev. of log of b
    real<lower=0> ri_kappa; // (Inv. sqrt) Clumping parameter in negative binomial
}

model {
    matrix[ntrap, nyear] mu; // Mean number of seeds per trap and year
    
    mu_disp ~ normal(p_mu_disp[1], p_mu_disp[2]);
    sd_disp ~ gamma(p_sd_disp[1], p_sd_disp[2]);
    mu_beta ~ normal(p_mu_beta[1], p_mu_beta[2]);
    sd_beta ~ normal(p_sd_beta[1], p_sd_beta[2]);
    
    beta_off ~ normal(0, 1);
    
    mu = calc_mu(mu_disp, sd_disp, beta_off, mu_beta, sd_beta, trap_area, nyear, 
                 ntree, ntrap, r, log_wgt_2pr, rmax, tree_size, size_density);
    
    nseed ~ poisson(to_vector(mu));
}

generated quantities {
    vector[ntrap*nyear] log_lik;
    int pval;
    int totseeds;
    int nnz;
    real rmode;
    real rmed;
    real rmean;
    
    {
        matrix[ntrap, nyear] mu;
        vector[ntrap*nyear] mu_v;
        int nseed_sim[ntrap*nyear];
        vector[ntrap*nyear] ll_sim;
        mu = calc_mu(mu_disp, sd_disp, beta_off, mu_beta, sd_beta, trap_area, nyear, 
                     ntree, ntrap, r, log_wgt_2pr, rmax, tree_size, size_density);
        mu_v = to_vector(mu);
        nnz = 0;
        for (i in 1:ntrap*nyear) {
            log_lik[i] = poisson_lpmf(nseed[i] | mu_v[i]);
            nseed_sim[i] = poisson_rng(mu_v[i]);
            ll_sim[i] = poisson_lpmf(nseed_sim[i] | mu_v[i]);
            if (nseed_sim[i] > 0) nnz = nnz + 1;
        }
        pval = sum(ll_sim) < sum(log_lik);
        totseeds = sum(nseed_sim);
    }
    
    rmode = exp(mu_disp - 2 * square(sd_disp));
    rmed = exp(mu_disp);
    rmean = exp(mu_disp + square(sd_disp) / 2);
}
