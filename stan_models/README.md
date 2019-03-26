# Stan models

Most Stan model files in this folder are named following the pattern: `disp_[disp_mod]_[err_mod]_[version].stan`, where *disp_mod* is the kernel used (2Dt, exppow, invpow, lognorm or weibull), *err_mod* is the count model used (nb or pois), and *version* is one of three model versions:

- `fixed`: Multiple years of seed data with same (fixed) tree distribution. Used only for the simulations (in `scripts_prior_checks` folder).

- `mat`: The count data is a matrix (traps x years). Used for seedling models where the same traps are active every year.

- `sparse`: The count data is a vector, with separate vector indicating the year index and trap index for each count (sparse matrix data). Used for seed and recruit models where the number of traps varies over study period.