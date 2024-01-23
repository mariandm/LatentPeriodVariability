# *Code for:* Accounting for Cellular-Level Variation in Lysis: Implications for Virus-Host Dynamics
*By Marian Dominguez-Mirazo, 2023*

## General Description

In this repository you can find all code needed to replicate the figures and analyses shown in the paper. Most code is written in Matlab R2021_a, except for the Parameter Inference framework, which is written in Julia v1.7.2. Information on the inference framework, including Julia versioning and dependencies, can be found in the README on folder *./ParameterInference. 

## Folder content description

### ./Data

Contains data retrieved from the literature and pre-run data to improve plotting speed. 

- **OneStep_data.csv :** Metadata table containing information on one-step growth curves retrieved from the literature. Required to plot Figure 1. A version of this table is available as Table S2. 

- **OneStep_datasets/ :** Contains the individual datasets for the one-step growth curves retrieved from the literature. Required to plot Figure 1. First column of the file is time, second column is free-virus measurement. The units for each column are shown in the metadata file. 

- **Simulated OneSteps :** One-step growth curve simulations for multiple coefficient of variation (CV) values for the three datasets used in parameter inference (see main text Table 2). The first column corresponds to the simulated CV, teh second column to the time of first visible burst. Required to plot Figure 5b. **Origin script: ./FigureScripts/GenerateOneStepTimesforFig5b.m**

- **host_CI_data1_6.mat :** Contains multiple total host simulations using the parameters inferred (see ParameterInference folder) for dataset: data1 and data ID: 6. Required to plot Figure 5a. **Origin plot: ./FigureScripts/GenerateCIforFig5a.m**

- **vir_CI_data1_6.mat :** Contains multiple free-virus simulations using the parameters inferred (see ParameterInference folder) for dataset: data1 and data ID: 6. Required to plot Figure 5a. **Origin plot: ./FigureScripts/GenerateCIforFig5a.m**

### ./FigureScripts

The folder contains all scripts required to replicate main and supplementary text figures, as well as some pre-run data required for plotting. 

- **Figure1_OneSteps.m**

- **Figure2_IndividualVariationEffect.m**

- **Figure3_identifiability.m**

- **Figure4_multicycle.m**

- **Figure5a_inferenceExample.m : ** The data required to plot this figure is obtained from the Parameter Inference framework. See the ParameterInference folder README for information on how to generate the data. 

- **Figure5b_inferenceAll.m :** The data required to plot this figure is obtained from the Parameter Inference framework. See the ParameterInference folder README for information on how the data was generated. 

- **FigureS1_EtoCV.m**

- **FigureS2_insilicodata.m :** The data required to plot this figure is obtained from the Parameter Inference framework. See the ParameterInference folder README for information on how the data was generated. 

- **FigureS3_otherParams.m : ** The data required to plot this figure is obtained from the Parameter Inference framework. See the ParameterInference folder README for information on how the data was generated. 

### ./Figures
Figure storage. Output folder for figure producing scripts found in ./FigureScripts folder. 

### ./ParameterInference
The folder contains all code required to run a parameter inference framework that predicts host and viral-life history traits based on parameter fitting. The code is written in Julia v1.7.2 and contains its own README. 


