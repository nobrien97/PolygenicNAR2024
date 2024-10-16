#!/bin/bash -l
#PBS -q hugemem
#PBS -l walltime=48:00:00
#PBS -l ncpus=4
#PBS -l mem=1470GB
# Purpose: Sets up a sql database for mutation data

# Path to updated version of sqlite
SQLITE3=~/Tools/sqlite/sqlite3
SAVEDIR=/storage_path/runSims/calcMutationStats
RUNDIR=/scratch_path

cd $PBS_O_WORKDIR

# Loads data into data frame and attaches indices
cat ./load_data.sql | $SQLITE3 $RUNDIR/epistasis.db 2>/dev/null
cat ./add_indices.sql | $SQLITE3 $RUNDIR/epistasis.db 2>/dev/null

mv $RUNDIR/epistasis.db* $SAVEDIR