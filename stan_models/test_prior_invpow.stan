functions {
    // theta: r, x_r: (a, k)
    vector find_med(vector y, vector theta, real[] x_r, int[] x_i) {
        vector[1] m;
        m[1] = pow(x_r[1], x_r[2] - 1) * (x_r[2] * y[1] + x_r[1]) /
            pow(y[1] + x_r[1], x_r[2]) - 0.5; 
        return m;
    }
}

data {
    real p_alpha[2];
    real p_k_real[2];
}

generated quantities {
    
    real alpha;
    real k_real;
    real k; 
    real a_k[2]; // Array containing a and k
    vector[1] guess;
    vector[0] theta;
    int x_i[0];
    real<lower=0> rmed;
    real rmean;

    alpha = normal_rng(p_alpha[1], p_alpha[2]);
    k_real = normal_rng(p_k_real[1], p_k_real[2]);
    k = inv_logit(k_real)*3 + 2;
    
    a_k[1] = exp(alpha);
    a_k[2] = k;
    guess[1] = a_k[1];
    
    rmed = algebra_solver(find_med, guess, theta, a_k, x_i)[1];
    rmean = 2 * a_k[1] / (a_k[2] - 2);
}
