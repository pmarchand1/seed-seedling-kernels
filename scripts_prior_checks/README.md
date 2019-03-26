# Prior distributions and model checks

Scripts in this folder serve to check model behavior with simulated data, prior to fitting the models to the actual study data.

- `prior_disp_stats.R`: Calculate the distribution of the median and mean dispersal distances corresponding to a given set of dispersal parameter prior distributions (using the `test_prior_[...].stan` files in the `stan_models` folder). This is useful to choose weakly informative priors that result in plausible dispersal statistics.

- `run_seed_model_sim_stan.R`: For a given species, dispersal model and count model, (1) draw a parameter set from the prior distributions, (2) simulate five years of seed dispersal based on trees present in 2010 census, (3) fit the model to the simulated data with Stan and (4) save the results in a file. The script is structured to be easily run in parallel on a computing cluster to simulate multiple runs of same model.

- `load_sim_data.R`: Called by script above to load the data for a simulation (tree census and seed trap locations).

- `analyze_sim_results.R`: Analyze the output from a set of simulations from `run_seed_model_sim_stan.R` to verify the consistency and precision of posterior estimates. This script generates the results found in Appendix B of the paper (*Prior distributions and model checks*). 