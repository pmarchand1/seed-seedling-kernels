// Seed dispersal model with inverse power kernel

functions {
    //  Calculate expected number of seeds per trap and year
    matrix calc_mu(real a, real k, row_vector beta_off, real mu_beta, real sd_beta,
                   data real trap_area, int nyear, int ntree, int ntrap, 
                   data matrix r, data matrix wgt, data vector rmax, 
                   data vector tree_size, data real size_density) {
        row_vector[nyear] b; // Fecundity parameter by year
        matrix[ntrap, ntree] disp_kern; // Dispersal kernel
        vector[ntrap] offplot; // Radial integral of disp. kern from rmax to inf.
        matrix[ntrap, nyear] mu;

        for (j in 1:ntree) {
            for (i in 1:ntrap) {
                disp_kern[i, j] = wgt[i, j] * (k-1)*k / (2*pi()*square(a)) * 
                                        pow(1 + r[i, j]/a, -1-k);
            }
        }
    
        for (i in 1:ntrap) {
            offplot[i] = a^(k-1) * (k*rmax[i] + a) * pow(rmax[i] + a, -k);    
        }
    
        b = exp(mu_beta + sd_beta * beta_off);
        
        mu = trap_area * (disp_kern * tree_size + offplot * size_density) * b;
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
    vector<lower=0>[ntree] tree_size;
    real<lower=0> size_density;
    
    // Number of seeds in traps
    int<lower=0> nseed[ntrap*nyear];
    
    // Hyperparameters (mean and std. dev.) for parameter priors
    real p_alpha[2];
    real p_k_real[2];
    real p_mu_beta[2];
    real p_sd_beta[2];
    real p_ri_theta[2];
}

parameters {
    real alpha; // (Log) of scale parameter
    real k_real; // Logit transform of shape parameter

    row_vector[nyear] beta_off; // Fecundity offset by year
    real mu_beta; // Mean log of b
    real<lower=0> sd_beta; // Std. dev. of log of b
    real<lower=0> ri_theta; // (Inv. sqrt) Clumping parameter in negative binomial
}

transformed parameters {
    real a; // Scale parameter
    real k; // Shape parameter
    real theta; // Clumping parameter in negative binomial
    a = exp(alpha);
    k = 3 * inv_logit(k_real) + 2;
    theta = inv_square(ri_theta);
}

model {
    matrix[ntrap, nyear] mu; // Mean number of seeds per trap and year
    
    alpha ~ normal(p_alpha[1], p_alpha[2]);
    k_real ~ normal(p_k_real[1], p_k_real[2]);
    mu_beta ~ normal(p_mu_beta[1], p_mu_beta[2]);
    sd_beta ~ normal(p_sd_beta[1], p_sd_beta[2]);
    ri_theta ~ normal(p_ri_theta[1], p_ri_theta[2]);
    
    beta_off ~ normal(0, 1);

    mu = calc_mu(a, k, beta_off, mu_beta, sd_beta, trap_area, nyear, ntree, 
                 ntrap, r, wgt, rmax, tree_size, size_density);
    
    nseed ~ neg_binomial_2(to_vector(mu), theta);
}

generated quantities {
    vector[ntrap*nyear] log_lik;
    int pval;
    int totseeds;
    int nnz;
    real rmode;
    real rmean;
    
    {
        matrix[ntrap, nyear] mu;
        vector[ntrap*nyear] mu_v;
        int nseed_sim[ntrap*nyear];
        vector[ntrap*nyear] ll_sim;
        mu = calc_mu(a, k, beta_off, mu_beta, sd_beta, trap_area, nyear, ntree, 
                     ntrap, r, wgt, rmax, tree_size, size_density);
        mu_v = to_vector(mu);
        nnz = 0;
        for (i in 1:ntrap*nyear) {
            log_lik[i] = neg_binomial_2_lpmf(nseed[i] | mu_v[i], theta);
            nseed_sim[i] = neg_binomial_2_rng(mu_v[i], theta);
            ll_sim[i] = neg_binomial_2_lpmf(nseed_sim[i] | mu_v[i], theta);
            if (nseed_sim[i] > 0) nnz = nnz + 1;
        }
        pval = sum(ll_sim) < sum(log_lik);
        totseeds = sum(nseed_sim);
    }
    
    rmode = 0;
    rmean = 2 * a / (k - 2);
}
