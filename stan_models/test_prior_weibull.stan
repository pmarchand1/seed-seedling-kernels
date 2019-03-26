data {
    real p_alpha[2];
    real p_inv_k[2];    
}

generated quantities {
    
    real alpha;
    real<lower=0> inv_k; 
    real a;
    real rmed;
    real rmean;
    
    alpha = normal_rng(p_alpha[1], p_alpha[2]);
    inv_k = gamma_rng(p_inv_k[1], p_inv_k[2]);
    a = exp(alpha - log(log2()) * inv_k);
    
    rmed = exp(alpha);
    rmean = a * tgamma(1 + inv_k);
}
