#!/bin/bash -l
#PBS -l walltime=5:00:00
#PBS -l ncpus=10800
#PBS -l mem=42750GB
#PBS -l jobfs=10000GB  
  
ECHO=/bin/echo
JOBNAME=runSims/calcLD

#
# These variables are assumed to be set:
#   NJOBS is the total number of jobs in a sequence of jobs (defaults to 1)
#   NJOB is the number of the current job in the sequence (defaults to 0)
#   For this job, NJOBS should = 1
  
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
    mkdir -p /scratch_path/${JOBNAME}R
    mkdir -p /storage_path/${JOBNAME}R
    mkdir -p $HOME/tests/$JOBNAME/done

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
SAVEDIR=/storage_path/${JOBNAME}R

export ncores_per_task=4
export ncores_per_numanode=12

# Calculate the range of parameter combinations we are exploring this job
# CAUTION: may error if CUR_TOT is not a multiple of PBS_NCPUS - untested
CMDS_PATH=$HOME/tests/$JOBNAME/PBS/cmdsR.txt
CUR_TOT=$(cat $CMDS_PATH | wc -l)
CUR_MIN=$(($NJOB*$PBS_NCPUS+1))
CUR_MAX=$((($NJOB+1)*$PBS_NCPUS))

if [ $CUR_MAX -gt $CUR_TOT ]; then
    CUR_MAX=$CUR_TOT
fi

sed -n -e "${CUR_MIN},${CUR_MAX}p" $CMDS_PATH > ./JOB_PATH.txt

# Run on your HPCs taskfarming software - mpirun, parallel etc.
cat ./JOB_PATH.txt | parallel -j $((PBS_NCPUS/ncores_per_task)) pbsdsh -n {%} -- bash -l -c '{}'

$ECHO "All jobs finished, moving output..."

# Combine output into a single file
cd /scratch_path/${JOBNAME}R/

cat ./out_LD_* >> $SAVEDIR/out_LD.csv
cat ./out_LDf_* >> $SAVEDIR/out_LDf.csv
cat ./out_LD_raw_* >> $SAVEDIR/out_LD_raw.csv

# Delete loose files with seed and model indices
find -regex ".*[0-9]*_*[0-9].csv+" -delete

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
    qsub -v NJOBS=$NJOBS,NJOB=$NJOB ./calcLDR.sh
else
    $ECHO "Finished last job in sequence of $NJOBS jobs"
fi