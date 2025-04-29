#!/bin/bash -l
#PBS -P ht96
#PBS -q normalsr
#PBS -l walltime=24:00:00
#PBS -l ncpus=520
#PBS -l mem=2500GB
#PBS -l jobfs=1000GB
#PBS -l storage=scratch/ht96+gdata/ht96

# Run with -v NJOBS=0

# Calculates statistics for mutation data
ECHO=/bin/echo
SUBJOBNAME=newMotifs/randomisedStarts
JOBNAME=epistasisDensity
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
    mkdir -p /scratch/ht96/nb9894/$TOTALJOBNAME
    mkdir -p /g/data/ht96/nb9894/$TOTALJOBNAME
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

# Pre-job file manipulation goes here ...
# 
# INSERT CODE
cd $PBS_O_WORKDIR
SAVEDIR=/g/data/ht96/nb9894/$TOTALJOBNAME


# ========================================================================
# .... USER INSERTION OF EXECUTABLE LINE HERE 
# ========================================================================
# Make sure we're at the right place so we can find the bash script to run

$ECHO "Running nci-parallel..."
# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load nci-parallel/1.0.0a
export ncores_per_task=1
export ncores_per_numanode=13           # 13 cores per NUMA node for sapphire rapids nodes

# Calculate the range of parameter combinations we are exploring this job
CMDS_PATH=$HOME/tests/$TOTALJOBNAME/PBS/cmds.txt
CUR_TOT=$(cat $CMDS_PATH | wc -l)
CUR_MIN=$(( ( ($NJOB*$PBS_NCPUS+1) + ($ncores_per_task-1) ) / $ncores_per_task ))
CUR_MAX=$(( ( ( ($NJOB+1) * $PBS_NCPUS ) + ($ncores_per_task-1) ) / $ncores_per_task ))

if [ $CUR_MAX -gt $CUR_TOT ]; then
    CUR_MAX=$CUR_TOT
fi

sed -n -e "${CUR_MIN},${CUR_MAX}p" $CMDS_PATH > ./JOB_PATH.txt

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file ./JOB_PATH.txt --timeout 172800

$ECHO "All jobs finished, moving output..."

# Combine output into a single file
cd /scratch/ht96/nb9894/$TOTALJOBNAME

cat ./d_epi_density* >> $SAVEDIR/d_epi_density.csv
cat ./d_epi_freqweight_density* >> $SAVEDIR/d_epi_freqweight_density.csv
cat ./d_epi_mean* >> $SAVEDIR/d_epi_mean.csv
cat ./d_epi_freqweight_mean* >> $SAVEDIR/d_epi_freqweight_mean.csv

# Without molecular components
cat ./d_epi_nomolcomp_density* >> $SAVEDIR/d_epi_nomolcomp_density.csv
cat ./d_epi_freqweight_nomolcomp_density* >> $SAVEDIR/d_epi_freqweight_nomolcomp_density.csv
cat ./d_epi_nomolcomp_mean* >> $SAVEDIR/d_epi_nomolcomp_mean.csv
cat ./d_epi_freqweight_nomolcomp_mean* >> $SAVEDIR/d_epi_freqweight_nomolcomp_mean.csv

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
    qsub -v NJOBS=$NJOBS,NJOB=$NJOB ./$JOBNAME.sh
else
    $ECHO "Finished last job in sequence of $NJOBS jobs"
fi 
