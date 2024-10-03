# Author: Nick O'Brien
# Runs analysis on data for the paper 
# The genetic architecture of polygenic adaptation under a network-derived trait
# O'Brien et al. 2024

# Replace these paths with the path to where you saved the 
# PolygenicNAR2024 repo and where you saved the dataset
repoPath <- "/path/to/PolygenicNAR2024"
dataPath <- "/path/to/data/"

setwd(repoPath)

# Load functions
source("./R/helperFunctionsAndSetup.R")

# Setup data
source("./R/wrangle_data.R")

# Analyse epistasis data
source("./R/mutationScreenExp.R")

# Analyse LD data
source("./R/figureFunctionsAndSetup.R")

# Analyse G matrices
source("./R/figures.R")
