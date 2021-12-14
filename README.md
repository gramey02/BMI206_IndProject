# BMI 206 Individual Project -- Last Updated 12/13/21
Repo for all files relating to my individual project for BMI 206, Fall 2021. Analysis was based on the paper that can be found here: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2897-z.

## Data files
Project was conducted on a directed network of E.coli metabolic reactions stored in an adjacency matrix (AdjData_ecoli.xlsx) and an E.coli metabolic network (iJO1366.mat) from the BiGG database which can be found here: http://bigg.ucsd.edu/models/iJO1366.

## .mat files
  - Ecoli_FBA.m - Script to run Flux balance analysis (FBA) to determine essential reactions
  - Ecoli_cascadeNums.m - Script to run the cascade algorithm on each node
  - cascade.m - algorithm to calculate cascade number for a node (from Kim et al., 2019)
  - in_degrees.m - algorithm to calculate in-degree for each node (from Kim et al., 2019)

## .py files
Please see the "compute graph metrics" file and the corresponding function scripts in "lib.py" of the code folder of the github site: https://github.com/eelxela/bmi206-project. These were scripts written during the group analysis project that I ran in order to calculate all graph centrality metrics (excluding cascade number).

## .R files
   - NewEcoliFig_Stats.R - Script to clean and merge output data from different sources (Python and MATLAB scripts) and create a figure showing 5 different centrality metrics and essential reaction proportions. This script was also used to determine if the correlation between cascade number and the % of essential reaction was significant which was my extension of the author's original figure.
   - CorCompare.R - Script to run a permutation test on the difference between correlation coefficients calculated with Pearson's method or Spearman's method. Final output is a plot showing how the observed difference for each centrality metric compares to a null distribution.
   - LogisticRegr.R - Script to construct a logistic regression model that can predict whether or not a reaction in the network is essential based on the non-binned values of that node's centrality metrics. Model performance metrics were also calculated in this script.
