# Kreg_GSA
 R version of K regulation model used for Morris Method analysis

 # Main driver files
These are the main driver files to compute simulation results and the Morris Method analysis 
**driver_singlemeals.r**
Runs a baseline simulation for each of the single meal experiments (Meal + KCl, Meal Only, KCl Only).

**runODEMorris_MealKCl.r**
Runs Morris method analysis for the Meal + KCl experiments.

**runODEMorris_MealOnly.r**
Runs Morris method analysis for Meal Only experiments

**runODEMorris_KClOnly.r**
Runs Morris method analysis for KCl Only experiments

**runMorris_SS.r**
Runs Morris method analysis for steady state results

 ## Results/
 Saved results from the simulations and global sensitivity analysis in this folder.
 Use **convert_MorrisResults.r** and **convert_MorrisSS.r** to convert .RData results files into .csv for plotting.

 ## Figures/
 Scripts for making figures from the results

 **manu_figs/** Copy of figures for manuscript

 **results_final/** simulation and Morris analysis results used for plotting

 **plot_MealSims.m** plot simulation results for the 3 single meal simulation types (Meal Only, KCl Only, Meal + KCl)

 **plot_MI_allmeal.m** plot Morris Indices for all 3 single meal simulation types

 **plot_MorrisSS.m** plot Morris analysis results for steady state

 **plot_Morris_KClOnly.m** plot Morris analysis results for KCl Only simulations

 **plot_Morris_MealOnly.m** plot Morris analysis results for Meal Only simulations

 **plot_Morris_MealKCl.m** plot Morris analysis results for Meal + KCl simulations