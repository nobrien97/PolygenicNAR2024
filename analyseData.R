# Author: Nick O'Brien
# Runs analysis on data for the paper 
# The genetic architecture of polygenic adaptation under a network-derived trait
# O'Brien et al. 2024

# Replace these paths with the path to where you saved the 
# PolygenicNAR2024 repo and where you saved the dataset
REPO_PATH <- "/path/to/PolygenicNAR2024/"
DATA_PATH <- "/path/to/data/"

setwd(REPO_PATH)

# Load functions
source("./R/helperFunctionsAndSetup.R")

# Setup data
source("./R/wrangle_data.R")

# Analyse epistasis data
source("./R/epistasis.R")

# Analyse LD data
source("./R/LD.R")

# Analyse G matrices
source("./R/Gmatrix.R")

# Find effect sizes of individual mutations
source("./R/FX.R")

# Phenotype figures
source("./R/adaptiveWalk.R")
