data {
    real p_alpha[2];
    real p_inv_k_real[2];
}

generated quantities {
    real alpha;
    real inv_k_real;
    real<lower=0, upper=2> inv_k; 
    real a;
    real k;
    real rmed;
    real rmean;
    
    alpha = normal_rng(p_alpha[1], p_alpha[2]);
    inv_k_real = normal_rng(p_inv_k_real[1], p_inv_k_real[2]);
    inv_k = inv_logit(inv_k_real) * 2;
    a = exp(alpha);
    k = 1/inv_k;
    
    rmed = sqrt(a * (pow(2, inv_k) - 1));
    rmean = sqrt(pi() * a) / 2 * tgamma(k - 0.5) / tgamma(k);
}
