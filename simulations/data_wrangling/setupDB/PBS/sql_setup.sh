#!/bin/bash -l
#PBS -l walltime=48:00:00
#PBS -l ncpus=48
#PBS -l mem=190GB

# Purpose: Sets up a sql database for mutation data

# Path to updated version of sqlite
SQLITE3=~/Tools/sqlite/sqlite3
SAVEDIR=/storage_path/runSims
RUNDIR=/scratch_path/

cd $PBS_O_WORKDIR

# Loads data into data frame and attaches indices
cat ./load_data.sql | $SQLITE3 $RUNDIR/standingVarMuts.db 2>/dev/null
cat ./add_indices.sql | $SQLITE3 $RUNDIR/standingVarMuts.db 2>/dev/null

mv $RUNDIR/standingVarMuts.db* $SAVEDIR