# SHAPPEing agent-based models to spatial data: Sequential Hooke’s law Algorithm for Point Pattern reconstruction with Environmental constraints
Authors: Blanche Mongeon and Morgan Craig

Corresponding author: blanche.mongeon@umontreal.ca

## Data sets from original article
### Tumour-microenvironments reconstructions
The folder `IMC_TMEReconstructions` contains:
1. the 87 original IMC samples,
2. their corresponding 10 same-size ECSR reconstructions (or CSR reconstructions for homogeneous samples),
3. their corresponding 10 same-size SHAPPE reconstructions,
4. their corresponding 10 _simplified_ SHAPPE reconstructions (i.e., without additional repulsion/attraction from exclusion/inclusion zones),
5. their corresponding 10 increased size SHAPPE reconstructions, and
6. results discussed in the article.

More specifically:

- File `SOCPatientsInfo.csv` contains relevant information for the 87 IMC samples (e.g., IMC size, number of cells of each cell type, etc.)
- File `inhomogeneityFinalResults_50.mat` contains the inhomogeneity identification results for the 87 samples.
  The first row refers to the first step (with the empty space function null enveloppe). The second row refers to the second step, where we apply the minimal area threshold.
- File `inhomogeneity_replicates.mat` contains the inhomogeneity status (1 if inhomogeneous, 0 if not) of the 87 samples (1st column) and of their 10 same-size SHAPPE reconstructions (columns 2 to 11).
- Files `RipleyLBs_markBootstrap.mat` and `RipleyUBs_markBootstrap.mat` contain respectively the lower and upper bounds of 95% confidence interval bounds of
  Ripley's H functions. First dimension is the sample, second and third are the cell types used in the evaluation of the function, and last is the distance (20, 40 and 60um).
- Files `insideBoundsRandom.mat`, `insideBoundsSHAPPE.mat` and `insideBoundsimplifiedSHAPPE.mat` contain the fraction of the 10
  same-size ECSR/CSR, SHAPPE and simplified SHAPPE that fall within the 95% confidence interval bounds of Ripley's H function at distances 20 (column 1), 40 (column 2) and 60um (column 3)
  for each IMC sample (indicated by the row). 
- File `hyperparametersResultTable_multipleScores_final.csv` contains the average normalised absolute clustering errors at distance 20, 40 and 60um for each hyperparameter combination.
- File `hyperparametersResultTable_combinedScores_final.csv` contains the mean of the average normalised absolute clustering errors at distance 20, 40 and 60um for each hyperparameter combination.
- File `hyperparameterImportance.csv` contains the normalised importance for each hyperparameters,
  with regards with the average normalised absolute clustering errors at distance 20, 40 and 60um, and their mean), as obtained with Optuna.
- Files `SHAPPE_HyperparameterOptimisation_SameSize.mat`, `SHAPPE_HyperparameterOptimisation_Simplified.mat` and `SHAPPE_HyperparameterOptimisation_DomainIncrease.mat`,
  and more specifically their `resultsOptimised` fields, contain the performance metrics (rows) for each of the 87 samples (columns). Metrics are in the order:
1. Excess out of the 95 C.I. of H(20um),
2. Excess out of the 95 C.I. of H(40um),
3. Excess out of the 95 C.I. of H(60um),
4. Raw difference with H(20um),
5. Raw difference with H(40um),
6. Raw difference with H(60um),
7. eDD,
8. eNN,
9. eSigma,
10. ePi.
  Fields `resultsCSR` and `resultsECSR` in the file `SHAPPE_HyperparameterOptimisation_SameSize.mat` contain the same metrics for the CSR and ECSR reconstructions, respectively.

The file `IMC_Results.m` allows any user to recreate the figures in the original article.

The code to generate the reconstructions is found in the main folder:

- `SHAPPE_commands.m` can be used to regenerate the SHAPPE reconstructions.
- `SHAPPE_commands_IMC_DomainEnlargement.m` contains the code to generate the larger domain reconstructions, and
  `SHAPPEPerformanceEvaluation_EnlargedDomains.m` contains the code to evaluate them.

File `SHAPPEResults.pdf` contains all the same-size SHAPPE reconstructions (with and without additional environmental constraints) and their corresponding metrics.

### PhysiCell

The ABM and IMC samples share the same cell types, but not the same cell type codes. Hence, the file `generateABMInitialConditions.m` map the larger domain 
reconstructions with IMC cell types code to aBM cell type codes using the `mapIMC_to_ABM.m` function. It also samples additional repulsion regions in an area 
equal to twice the larger domain area, to account for domain increase during the ABM simulations.

The folder `PhysiCell_SHAPPE` contains the requirements, the code and the scripts for the ABM simulations.

Subfolder `scripts` contains the scripts used on Digital Research Alliance of Canada servers. 

## Using SHAPPE

SHAPPE can be used via the file `SHAPPE_commands.m`. A user-guide is currently under preparation to facilitate its use and possible reimplementations.

## Additional information

Any additional information and/or data requests can be sent to blanche.mongeon@umontreal.ca .
