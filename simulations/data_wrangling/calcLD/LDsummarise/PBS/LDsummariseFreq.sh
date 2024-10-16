#!/bin/bash -l
#PBS -l walltime=5:00:00
#PBS -l ncpus=10800
#PBS -l mem=42750GB
#PBS -l jobfs=4500GB
# Run with -v NJOBS=0

# Calculates statistics for mutation data
ECHO=/bin/echo
SUBJOBNAME=runSims
JOBNAME=calcLD/LDsummarise
TOTALJOBNAME=$SUBJOBNAME/$JOBNAME
#
# These variables are assumed to be set:
#   NJOBS is the total number of jobs in a sequence of jobs (defaults to 1)
#   NJOB is the number of the current job in the sequence (defaults to 0)
  
if [ X$NJOBS == X ]; then
    $ECHO "NJOBS (total number of jobs in sequence) is not set - defaulting to 1"
    export NJOBS=1
fi
  
if [ X$NJOB == X ]; then
    $ECHO "NJOB (current job number in sequence) is not set - defaulting to 0"
    export NJOB=0
    # Since this is the first iteration, create our folders
    $ECHO "Creating outputs folders..."
    cd $PBS_O_WORKDIR

    # Make output folder
    mkdir -p /scratch_path/$TOTALJOBNAME
    mkdir -p /storage_path/$TOTALJOBNAME
    mkdir -p $HOME/tests/$TOTALJOBNAME/done

fi

#
# Quick termination of job sequence - look for a specific file 
#
if [ -f STOP_SEQUENCE ] ; then
    $ECHO  "Terminating sequence at job number $NJOB"
    exit 0
fi

$ECHO "Starting job $NJOB of $NJOBS"

cd $PBS_O_WORKDIR
SAVEDIR=/storage_path/$TOTALJOBNAME


export ncores_per_task=1
export ncores_per_numanode=12         

# Calculate the range of parameter combinations we are exploring this job
CMDS_PATH=$HOME/tests/$TOTALJOBNAME/PBS/cmds_freq.txt

# Run on your HPCs taskfarming software - mpirun, parallel etc.
cat ./JOB_PATH.txt | parallel -j $((PBS_NCPUS/ncores_per_task)) pbsdsh -n {%} -- bash -l -c '{}'

$ECHO "All jobs finished, moving output..."

# Combine output into a single file
cd /scratch_path/$TOTALJOBNAME

cat ./sumLD* >> $SAVEDIR/sumLD_df_new.csv

# Delete loose files with seed and model indices
find -regex ".*[0-9]+.csv" -delete

# 
# Check the exit status
#
errstat=$?
if [ $errstat -ne 0 ]; then
    # A brief nap so PBS kills us in normal termination
    # If execution line above exceeded some limit we want PBS
    # to kill us hard
    sleep 5 
    $ECHO "Job number $NJOB returned an error status $errstat - stopping job sequence."
    exit $errstat
fi

#   
# Are we in an incomplete job sequence - more jobs to run ?
#   
if [ $NJOB -lt $NJOBS ]; then
# Now increment counter and submit the next job
# 
    NJOB=$(($NJOB+1))
    $ECHO "Submitting job number $NJOB in sequence of $NJOBS jobs"
    cd $PBS_O_WORKDIR
    qsub -v NJOBS=$NJOBS,NJOB=$NJOB ./LDsummariseFreq.sh
else
    $ECHO "Finished last job in sequence of $NJOBS jobs"
fi 
