# PolygenicNAR2024
Simulation and analysis scripts for O'Brien et al. 2024 - The genetic architecture of polygenic adaptation under a network-derived trait

## Running analysis
The analysis can be done on a reasonably modest computer, but you will need a fair amount of RAM to complete the mutational effect analysis. Somewhere in the range of 16GB should be adequate. 
You can download the data at (link). After downloading the data, open analyseData.R and edit REPO_PATH and DATA_PATH to match where you have saved this repo and the data. You can then run analyseData.R to reproduce the figures from the paper. 

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

