functions {
    // theta: r, x_r: (a, k)
    vector find_med(vector y, vector theta, real[] x_r, int[] x_i) {
        vector[1] m;
        m[1] = gamma_cdf(pow(y[1] / x_r[1], x_r[2]), 2/x_r[2], 1) - 0.5; 
        return to_vector(m);
    }
}

data {
    real p_alpha[2];
    real p_inv_k[2];    
}

generated quantities {
    
    real alpha;
    real<lower=0> inv_k; 
    real a_k[2]; // Array containing a and k
    vector[1] guess;
    vector[1] theta;
    int x_i[1];
    real rmed;
    real rmean;

    alpha = normal_rng(p_alpha[1], p_alpha[2]);
    inv_k = gamma_rng(p_inv_k[1], p_inv_k[2]);
    a_k[1] = exp(alpha - inv_k);
    a_k[2] = 1/inv_k;
    guess[1] = a_k[1];
    theta[1] = 0;
    x_i[1] = 0;
    
    rmed = algebra_solver(find_med, guess, theta, a_k, x_i)[1];
    rmean = a_k[1] * tgamma(3 * inv_k) / tgamma(2 * inv_k);
}
