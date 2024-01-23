# Parameter Inference and Data Generation Framework
*By Marian Dominguez-Mirazo, 2023*

## General description

This folder contains instructions on how to run a parameter inference framework that predicts host and viral life-history traits based on fitting 'multi-cycle response curves' to find the parameter set that reduces the error between simulations and data, assuming 3 replicates with noise. 
For details, please visit the Methods section and supplementary material.
This folder contains all the pre-run outputs necessary to replicate figures 5, S2 and S3. 

## Version and dependencies
This code was written in Julia v1.7.2

Packages required:

- CSV v0.10.4
 
- DataFrames v1.3.4

- Dierckx v0.5.2

- DifferentialEquations v7.2.0

- Distributed 

- Distributions v0.25.67

- DSP v0.7.7

- ForwardDiff v0.10.35

- GLM v1.8.0

- JLD v0.13.2

- MCMCChains v5.3.1

- ParallelDataTransfer v0.5.0

- QuadGK v2.5.0

- SparseGrids v2.0.1

- SpecialFunctions v2.1.7

- StatsBase v0.33.21

- Turing v0.21.10






## Folder content description

### ./CreateData
Creates 3 replicates from sampled simulations with added noise for 3 parameter sets (see Methods section, Table 2), and its corresponding no-phage control for host-only parameter estimation. 
Each parameter set is ran for 10 coefficient of variation (CV) values, resulting in 30 simulations with 3 replicates each. 

#### Data IDs

|Parameter set | Set ID |
| :-----  | :-------|
|E.coli and lambda	  | data1 |
|P. marinus and PHM2  | data2 |
|E.hux and EhV			  | data3 |

#### Simulation IDs

|CV value 	|Simulation ID|
|:--- | :---|
|0.5 		|1|
|0.45		|2|
|0.4 		|3|
|0.35		|4|
|0.3 		|5|
|0.25 	|6|
|0.2 		|7|
|0.15 	|8|
|0.1 		|9|
|0.05 	|10|

*Script:*

- **createData.jl** : Generates all the following output files

*Output:*

- **{dataID}/novir_noise_{simulationID}** : no-phage control replicates, first column are timepoints (hr), other columns are total host replicates (CFU/ml)

- **{dataID}/novirOG_{simulationID}** : no-phage control, original simulation, first column are timepoints (hr), second column is Susceptible cells

- **{dataID}/vir_noiseh_{simulationID}** : multi-cycle curve replicates for total host counts, first column are timepoints (hr), other columns are total host replicates (CFU/ml)

- **{dataID}/vir_noisev_{simulationID}** : multi-cycle curve replicates for viral counts, first column are timepoints (hr), other columns are viral replicates (PFU/ml)

- **{dataID}/virOG_{simulationID}** : multi-cycle curve, original simulation, first column are timepoints (hr), second column is Susceptible cell counts (CFU/ml), intermediate columns are intermediate infected states (CFU/ml), last column are viral counts (PFU/ml)

- **dataID/parsOG** : Parameters with which data was generated, order is mu, K, phi, burst size, eta, CV.

- **{dataID}/x0OG** : Initial conditions with which data was generated, first value are susceptible cells (CFU/ml), second value is free virus (PFU/ml)

### ./HostOnlyParams

*Script: *

- **host_inference.jl** : Infers host-only parameters (mu, growth rate and K, carrying capacity), using host-only growth replicates from the CreateData folder in 2 steps. 

	Step 1) Feature extraction from growth curve
	
	Step 2) MCMC with priors informed by step 1)
	
*Output:*

- **{dataID}/prior_estimates_{simulationID}** : Initial parameter estimations calculated from the growth slope and plateau intersect 

- **{dataID}/hostchain_{simulationID}** : MCMC chains, first column corresponds to mu, second column to K, third column to the error chain.

### ./VirusHostParams

*Script:*

- **step1_likelihood.jl** : Scripts that searches parameter spaces to get the Maximum Likelihood parameter combination

*Output: *

- **step1_likelihood/{dataID}/estimates_{simulationID}** : Maximum Likelihood parameter combination in the order mu (imported from HostOnlyParams output), K (imported from HostOnlyParams output), phi, beta, eta, CV

- **step1_likelihood/{dataID}/ll_{simulationID}** : Likelihood values associated to explored parameter combinations, first column corresponds to eta, second column to CV, thrid column to phi, fourth column to beta, and fifth column to the likelihood value

*Script:*

- **step2_round1_MCMC.jl** : First round of MCMC for virus-host parameter inference, prior distributions are informed from the step1_likelihood output

*Output:*

- **step2_MCMC/round1/{dataID}/viruschain_{simulationID}_round1** : Chains for first round of MCMC, first column corresponds to log10(phi), second column to beta, third column to eta, fourth column to CV, fifth and sixth to host and virus likelihood chains. The script was originally run in a computer cluster, individual chains were manually concatenated.

*Script:*

- **step2_round2_MCMC.jl** : Second round of MCMC for virus-host parameter inference, prior distributions are informed from the step2_round1_MCMC output

*Output:*

- **step2_MCMC/round1/{dataID}/viruschain_{simulationID}_round2** : Chains for second round of MCMC, first column corresponds to log10(phi), second column to beta, third column to eta, fourth column to CV, fifth and sixth to host and virus likelihood chains. The script was originally run in a computer cluster, individual chains were manually concatenated.

### ./Functions

This folder contains functions used across main scripts

- **ODE_SEnIV.jl** : ODE function for julia solver, see Methods for model description
- run_vir.jl : Function to run multi-cycle response curve simulations. It approximates CV to the closest value described by an integer number of E compartments (E), see Figure S1 for details.

- **segment_estimate.jl** : Step 1 of host-only parameter estimation. Gets point estimates for growth rate and carrying capacity using host-only growth curves.

- **prior_estimate_ll.jl** : Main function for Step 1 of virus-host parameter inference

- **likelihood_matrix.jl** : Runs simulations of multiple parameter value combinations

- **get_likelihood.jl** : Calculates likelihood for a matrix of simulations against data replicates


## Example run

```
dataid="data1"
id="1"

# Create data
cd ./createData
julia ./CreateData.jl dataid 

# Estimate host-only parameters
cd ../HostOnlyParams
julia HostOnlyParams/host_inference.jl dataid id

# Step 1 to estimate virus-host parameters
cd ../VirusHostParams
julia VirusHostParams/step1_likelihood.jl dataid id

# Step 2 to estimate virus-host parameters

# First round of MCMC
julia VirusHostParams/step2_MCMC_round1.jl dataid id

# Second round of MCMC
julia VirusHostParams/step2_MCMC_round2.jl dataid id

```

The final parameter estimates will be contained in the chains of the host-only parameter inference (file HostOnlyParams/dataID/hostchain_simulationID, columns 1 and 2), and the second round of MCMC (file VirusHostParams/step2_round2_MCMC/dataID/viruschain_simulationID, columns 1 to 4).

This code was originally ran on a computer cluster. 
