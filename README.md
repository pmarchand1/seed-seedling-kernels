# Code for estimation of seed, recruit and seedling dispersal models using data from the Barro Colorado Island (BCI) forest dynamics plot


## Terminology and abbreviations

- The study fits dispersal models to data from three life stages ("stage"): seeds, new seedling recruits ("recruits") and 20-cm seedlings ("seedlings").

- The dispersal kernels include the 2D *t* ("2Dt"), exponential power ("exppow"), inverse power ("invpow"), log-normal ("lognorm") and Weibull ("weibull") distributions.

- The count models (relative to expected counts per trap) include the Poisson ("pois") and negative binomial ("nb") distributions.

- Since the same Stan models are used for each life stage, "seed" in the model can designate either seeds, recruits or seedlings, and "trap" can designate either seed traps or the sampling quadrats for recruits and seedlings.


## Content overview

Main folder files:

- `run_single_model.R`: For a given species, life stage, dispersal model and count model, (1) create the appropriate data subset, (2) fit the appropriate hierarchical Bayesian model to the data with Stan and (3) save the results in a file. The script is structured to be easily run in parallel on a computing cluster, with command line arguments to indicate which species/stage/model to run.

- `calc_dists_weights_funcs.R`: Functions to calculate tree-trap distances and edge-correction weights to account for trees outside plot boundaries.

- `data_subset_funcs.R`: Functions that process the data files to get trees, seeds, recruits or seedlings for a given species and year, to serve as input to the models.

- `disp_kernel_funcs.R`: Functions to calculate dispersal kernel-related quantities (probability density, cumulative probability, summary statistics) for each model.

- `disp_priors.csv` and `repr_priors.csv`: Fixed parameters used to define the prior distributions for the base parameters of the hierarchical models.

Subfolders:

- `scripts_prior_checks`: Scripts to visualize prior distributions and perform simulation-based model checks.

- `scripts_post_process`: Scripts to analyze model outputs.

- `stan_models`: All the .stan model files.

Each folder includes its own README file.


## Structure of the input data

This section describes the data files used to run the scripts (not included in this repository). Except for the species attributes table (`sp_tab.csv`) in the main folder, all data files were located in a `data`subfolder. In the descriptions below, column names in each file are indicated in italics.

- `tree_census_intervals.csv`: Generated from the [BCI Forest Census Plot Data](https://repository.si.edu/handle/10088/20925), this dataset contains four fields identifying a tree (*treeID*, 6-letter lowercase species code *sp*, coordinates *x* and *y*), and six fields giving the *census_year*, *status* (A = alive, D = dead) and basal area (*ba*, in cm^2) for two consecutive censuses where the tree was recorded (*census_year1*, *census_year2*, etc.). Therefore, there is one row for each inter-census interval for each tree. This particular structure facilitates interpolation of basal area and survival between census years.

- `seeds_merged_calyr.csv`: Generated from seed surveys data at BCI (managed by S.J. Wright), this dataset contains columns for the trap ID (*trap*), trap coordinates *X* and *Y*, 4-letter uppercase species code (*sp*), the *year* and number of seed equivalents (*seed_eq*). There is one row for each species, trap and year (zero counts included). Not all traps were surveyed each year.

- `seed_traps.csv`: For convenience, this file records the position of each seed trap (columns *TRAP*, *X* and *Y*).

- `recruits_merged.csv`: Generated from recruit surveys data at BCI (managed by S.J. Wright), this dataset contains columns for the ID of the recruit survey *quadrat* and its coordinates (*X*, *Y*), the ID of the associated seed *trap*, the 4-letter uppercase species code (*sp*), the seed dispersal year (*seedyr*) and the number of new recruits in the quadrat (*recruit*). There is one row for each species, quadrat and year (zero counts included). Not all quadrats were surveyed each year.

- `recruit_traps.csv`: For convenience, this file records the position of each recruit quadrat (columns *trap*, *quadrat*, *X* and *Y*).

- `seedlings_tidy.csv`: Generated from 20-cm and taller seedling surveys at BCI (managed by L.S. Comita), this dataset contains one row by seedling and year, with columns for the seedling ID (*TAGF*), the 1-m^2 quadrat coordinates (*PX*, *PY*) where it was observed, the 6-letter uppercase species code (*SPP*), the *year* and *status* (A = alive, P = not yet observed).

- `sp_tab.csv`: This table has one row by species and contains the species name (*sp_name*), its 4-letter and 6-letter codes used in the BCI data (*sp_code4*, *sp_code6*), the estimated minimum basal area in cm^2 for reproductive trees (*rba_cm2*), and the maximum number of reproductive trees of that species at any time in the study period (*tree_count*). It also includes total counts of seeds, recruits and seedlings of that species in the data subset used for this study (*seed_count*, *recruit_count* and *seedling_near_traps*), as well as the number of trap-years or quadrat-years with non-zero counts (*seed_trap_yrs*, *recruit_trap_yrs* and *seedling_near_traps_yrs*). Species names and codes are based on the [BCI taxonomy file](https://repository.si.edu/handle/10088/32990).
