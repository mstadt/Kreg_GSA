# Kreg_GSA
 R version of K regulation model used for Morris Method analysis

 ## Results/
 Save results from the simulations and global sensitivity analysis in this folder.
 Use **convert_MorrisResults.r** and **convert_MorrisSS.r** to convert .RData results files into .csv for plotting.

 ## Figures/
 Scripts for making figures from the results
 **results_final** copy of results used to make the figures

 **plot_MealSims.m** plot simulation results for the 3 single meal simulation types (Meal Only, KCl Only, Meal + KCl)

# Main drivers
These are the main driver files
## driver_singlemeals.r
Runs a baseline simulation for each of the single meal experiments ("Meal + KCl", "Meal Only", "KCl Only").

## runODEMorris_MealKCl.r
Runs Morris method analysis for the "Meal + KCl" experiments.

## runODEMorris_MealOnly.r
Runs Morris method analysis for "Meal Only" experiments

## runODEMorris_KClOnly.r
Runs Morris method analysis for "KCl Only" experiments

## runMorris_SS.r
Runs Morris method analysis for steady state results