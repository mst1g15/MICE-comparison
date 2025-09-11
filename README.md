# ML-MICE comparison
Welcome! This repository contains R scripts to run simulation studies to compare the default and ML implementations of MICE. 

## Code for simulation - Single Covariate and Interaction settings

The repository contains the following R scripts: 

- `00_init.R` loads packages and settings.

- `01_generate_complete_data.R` contains a function for the data generating mechanism

- `01_simulate_missing_data.R` contains a function which induces missingness

- `02_impute_functions.R` contains a function which handles missingness with MICE using default and ML approaches 

- `03_simulation_function.R` contains the main wrapper function to run the simulation study 

- `04_run_sim.R` demonstrates how to run the simulation on the HPC 

## Code for simulation study based on the INOPulse trial 

- `07_mvpa_functions.R` contains functions to run the simulation study
  
- `07_mvpa_simulation_function.R` contains the main wrapper function to run the simulation study 

- `07_run_mvpa_sim.R` demonstrates how to run the simulation on the HPC 

## Results

Plots displaying results from the simulation studies are provided in the results folders: 

- Results_summarised contains plots presented in the main text of the paper and Supplementary File 2.

- Results_Single_Covariate and Results_Interactions contains plots with detailed results. 


Please contact Mia.Tackney@mrc-bsu.cam.ac.uk for any queries/comments.
