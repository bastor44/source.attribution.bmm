# source.attribution.bmm
## Project Description 
This repository was created in conjunction with a research project completed as part of the MSc in Statistics at Imperial College London in June-August 2025. This project aimed to adapt an existing Bayesian mixture model (BMM) for the identification of infectious disease transmission pairs by 1) evaluating the accuracy of the existing model, and 2) exploring alternative distributions for the  This project used Blenkinsop, et al. (2025) as a basis for both the models developed and the code included here. 

### Folder Organisation
- **data_Ams**: contains anonymised Amsterdam pairs data used in the study as an illustration of applying the proposed model to real-world data
- **data_other**: contains the simulated network of transmission events used in the simulation studies of the paper. Subsamples of the network are used to meet the mixing parameters set by the user
- **molecular clock**: figures of the molecular evolutionary clock model fit to the known Belgian transmission chain data (see data availability statement in associated paper)
- **R**: contains scripts with functions used in the analysis
    - formulate-Amsterdam-pairs.R - formulate phylogenetically possible pairs using epidemiological and clinical data. This script cannot be run in full due to data confidentiality.
    - molecular_clock_model_Gamma_hier.R - fit evolutionary clock model to known transmission chain data 
    - functions_plotting.R 
    - functions_simulation_scenarios.R
    - component_EM.R - function for running component-wise Expectation Maximisation (as laid out in Figueiredo & Jain, 2002)
    - test-CEM.R
- **scripts**: contains scripts to run the analysis on simulated data or observed Amsterdam data
    - make-stan-data.R
    - make-stan-data-Amsterdam.R
    - post-processing.R - produces convergence diagnostics, summarises generated quantities, including transmission flows for model fitted to simulated data
    - post-processing-Amsterdam.R - produces convergence diagnostics, summarises generated quantities, including transmission flows for model fitted to Amsterdam data
    - 2-stage-process.R
    - make-plot-simulated-networks.R
    - make-panel-plots-new.R
- **stan_model_files**:
    - mm_bgUnif_piGP_221027b.stan
    - mm_bgGMM_piGP_240711.stan
    - clock_model_gamma_hier_220315.stan
- **figures**

When creating simulated data, a new folder called **out_simulated** will be created, and folders within it will store the outputs of the simulations. When running models applied to the Amsterdam pairs data, the folder **out_Amsterdam** will be created, and the new output files will be stored within it. 

### Running Code 
To recreate the analyses done for this project, the arguments for input and output directories will need to be modified for the user. 

1. molecular_clock_model_gamma_hier_220315.stan to obtain the fitted values of the evolutionary clock model.
2. formulate-Amsterdam-pairs.R - run file setup, then lines 526-571 (save molecular clock quantiles)
3. make-stan-data.R (for Uniform background or GMM background fit to all pairs), 2-stage-process.R (for GMM background fit after background subtraction), or make-stan-data-Amsterdam.R, depending on whether applying the model to simulated data or observed data. Also, need to modify the corresponding stan model file to comment out or activate the input data and generated that requires knowledge of true transmission status (idx_true_pairs, true_flows, AE_from, and MAE_from)
4. post-processing.R or post-processing-Amsterdam.R and make-panel-plots-new.R

### Sources 
- Blenkinsop, A., Sofocleous, L., Lauro, F. D., Kostaki, E. G., Bezemer, D., Reiss, P., Pantazis, N., & Ratmann, O. (2025). Bayesian mixture models for phylogenetic source attribution from consensus sequences and time since infection estimates. Statistical Methods in Medical Research. https://doi.org/10.1177/09622802241309750.
- M. A. T. Figueiredo and A. K. Jain, "Unsupervised learning of finite mixture models," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 24, no. 3, pp. 381-396, March 2002, doi: 10.1109/34.990138.

