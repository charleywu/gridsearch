# Gridsearch
Complete code and data for replicating the results from:

Wu, C. M., Schulz, E., Speekenbrink, M., Nelson, J. D., & Meder, B. (2018). Generalization guides human exploration in vast decision spaces. Nature Human Behaviour. https://doi.org/10.1038/s41562-018-0467-4

An earlier version of this project (Experiment 2D and early modeling results) was released as:

Wu, C.M., Schulz, E., Speekenbrink, M., Nelson, J.D., Meder, B. (2017). Mapping the unknown: The spatially correlated multi-armed bandit. In G. Gunzelmann, A. Howes, T. Tenbrink, & E. J. Davelaar (Eds.), Proceedings of the 39th Annual Conference of the Cognitive Science Society (pp. 1357-1362). Austin, TX: Cognitive Science Society. doi: https://doi.org/10.1101/106286


## Getting Started and prerequisites

The environments used in Experiments 1 and 2 were generated using `samplePrior1D.py` and `samplePrior2D.py`, which requires numpy, GPy, and Scikit-learn. Experiment 3 uses a variety of agricultural data sampled from the R package `agridat` (see SI for full details), which are saved in `/experiment3/environments/agridat.json`.

99.9% of the rest of the data analysis uses R, with code separated into `\analysis1D`, `\analysis2D`, and `analysis3`. Required packages are specified at the top of each file. Model recovery code has it's own folder called \modelRecovery. Additionally, mismatch simulations are in \mismatch, where empirical simulations for each experiment are in R, while the generalized mismatch simulation is in python using the bayesian optimization library (https://github.com/jmetzen/bayesian_optimization). Mismatch simulations for Experiment 3 are in `\analysis3\mismatch.R`. 


## Authors

See the list of [contributors](https://github.com/charleywu/gridsearch/blob/master/contributors.txt) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
