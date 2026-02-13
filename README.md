# PolygenicNAR2024
Simulation and analysis scripts for O'Brien et al. 2024 - The genetic architecture of polygenic adaptation under a network-derived trait

## Running analysis
The analysis can be done on a reasonably modest computer, but you will need a fair amount of RAM to complete the mutational effect analysis. Somewhere in the range of 16GB should be adequate. 
You can download the data [here](https://espace.library.uq.edu.au/view/UQ:1a19d80). After downloading the data, open analyseData.R and edit REPO_PATH and DATA_PATH to match where you have saved this repo and the data. You can then run analyseData.R to reproduce the figures from the paper. 

## Running simulations
To rerun the simulations and produce the data linked above from scratch you will need access to an HPC: there are 21600 simulations to run, each takes 8-13hrs CPU time. The analysis is also quite expensive, especially the epistasis and LD calculations which have numerous steps.
In the PBS scripts, you will have to change the paths/setup details to suit your HPC environment. If your HPC is Slurm based, you will have to convert the PBS scripts over as well.

The order of the scripts/jobs you need to submit are:

- runSims (runs the SLiM simulations)
- setupDB (puts trait data and mutational data into a SQLite database)
- calcMutationStats (calculates pairwise epistasis, effect sizes, SFS, change in phenotype over time)
    - epistasis_sql_setup (puts epistasis data into a SQLite database)
- epistasisDensity (simplifies epistasis data into probability densities, means and standard deviations)
- getVA (calculates VA in Z expression and between molecular components)
- calcLD (reruns simulations to capture allele frequencies)
  - calcLDR (calculates LD using the above allele frequencies)
  - LDsummarise (average LD across groups)

## Data

This dataset contains files from O'Brien et al. "The genetic architecture of polygenic adaptation

under a network-derived trait". Below is a list of the data files and a description of what they summarise:

Modelindex is the row in combos.csv (i.e. which model is being run)

- standingVar_seeds.csv: set of 48 randomly sampled seeds used in the SLiM simulations.

- combos.csv: All 450 models ran in the SLiM simulations. Column 1 is the number of loci, column 2 is the mutational effect size variance, column 3 is the genome-wide recombination rate (per-locus/per-generation), column 4 is the model type (ODE is the K- model, K is K+).

- slim_qg.csv: trait data measured during the SLiM simulations. Columns: "generation", "seed", "modelindex", "mean heterozygosity", "additive variance", "mean phenotype", "phenotypic variance", "distance to the optimum", "mean fitness", "change in mean phenotype over the last 50 generations", "change in fitness over the last 50 generations", "mean aZ", "mean bZ", "mean KZ", "mean KXZ".

- d_fx.csv: Individual fitness effects of mutations. Measured in heterozygous state. Columns: "generation", "seed", "modelindex", "mutation type (e.g. aZ, bZ)", "mutation ID", "s (fitness effect)".

- d_epi_mean.csv: Mean and standard deviation per model of pairwise trait and fitness epistasis. Columns: "Progress to the optimum (0 - 100%)", "modelindex", "mean trait epistasis", "trait epistasis standard deviation", "mean fitness epistasis", "fitness epistasis standard deviation", "number of pairs evaluated for this model".

- d_epi_mean_percomp.csv: Mean and standard deviations of pairwise trait and fitness epistasis measured per molecular component. Columns: "Progress to the optimum (0 - 100%)", "modelindex", "molecular component", "mean trait epistasis", "trait epistasis standard deviation", "mean fitness epistasis", "fitness epistasis standard deviation", "number of pairs evaluated for this model/molecular component combination".

- d_epi_freqweight_mean_percomp.csv: Mean and standard deviations of pairwise trait and fitness epistasis measured pre molecular component, with pairs of loci having minor allele frequencies within 10% of each other. Columns are the same as d_epi_mean_percomp.csv.

- out_LD.csv: Pairwise LD summary per simulation. Columns: "generation", "seed", "modelindex", "mean LD (D)", "D standard deviation (D)", "total samples", "total samples where D > 0", "total samples where D < 0", "total samples where D ~ 0", 21 bins counting samples where D is in a range of 0.05 (e.g. 0 to 0.05, 0.2 to 0.25)

- out_LDf.csv: Pairwise LD summary per simulation, with pairwise comparisons between loci where the minor allele frequencies of both are within 10% of each other. Columns are the same as out_LD.csv

- d_VA_Z.csv: Additive variance/covariance estimates for molecular components/Z. Columns: "generation","seed","modelindex","VA_Z","VA_a","VA_b","VA_KZ","VA_KXZ","CVA_Z_a","CVA_Z_b","CVA_a_b","CVA_Z_KZ","CVA_a_KZ","CVA_b_KZ","CVA_Z_KXZ","CVA_a_KXZ","CVA_b_KXZ","CVA_KZ_KXZ","h2_Z","h2_a","h2_b","h2_KZ","h2_KXZ","method","optPerc","isAdapted","model","nloci","tau","r". VA is additive variance, CVA is additive covariance between two given components. h2 is narrow-sense heritability, method is the kernel method used (mkr or mrr from the bWGR R package), model, nloci, tau, r are model parameters from combos.csv.

- d_VA_noZ.csv: Additive variance/covariance estimates for molecular components with Z excluded from analysis. Columns are the same as d_VA_Z.csv, minus columns describing VA/CVA in Z.

- d_bootPCASim.RDS: R object containing the results from the PCA similarity factor experiment comparing G matrices.

- betareg_pcaSim_big.RDS: R object containing a beta regression comparing model types on PCA similarity factor.

- betareg_benmut.RDS: R object containing a beta regression comparing molecular components and models in the proportion of beneficial mutations.

- betareg_benmut_allmuttypes.RDS: R object containing a beta regression comparing only models in the proportion of beneficial mutations.

