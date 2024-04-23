#!/bin/bash -l
#PBS -P ht96
#PBS -l walltime=5:00:00
#PBS -l ncpus=9840
#PBS -l mem=19680GB
#PBS -l jobfs=19680GB
#PBS -l storage=scratch/ht96+gdata/ht96

SUBJOBNAME=standingVar
JOBNAME=getH2
TOTALJOBNAME=$SUBJOBNAME/$JOBNAME
TASKS_THIS_ITERATION=5
ECHO=/bin/echo

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
    mkdir /scratch/ht96/nb9894/${TOTALJOBNAME}
    mkdir /g/data/ht96/nb9894/${TOTALJOBNAME}
    mkdir $HOME/tests/${TOTALJOBNAME}/done

fi

# Quick terminate job sequence
if [ -f STOP_SEQUENCE ] ; then
    $ECHO  "Terminating sequence at job number $NJOB"
    exit 0
fi

cd $PBS_O_WORKDIR
SAVEDIR=/g/data/ht96/nb9894/${TOTALJOBNAME}

# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load nci-parallel/1.0.0a
export ncores_per_task=1
export ncores_per_numanode=12

# Calculate the range of parameter combinations we are exploring this job

# CAUTION: may error if CUR_TOT is not a multiple of PBS_NCPUS - untested
CMDS_PATH=$HOME/tests/${TOTALJOBNAME}/PBS/cmds.txt
CUR_TOT=$(cat $CMDS_PATH | wc -l)
CUR_MIN=$(($NJOB*$PBS_NCPUS*$TASKS_THIS_ITERATION+1))
CUR_MAX=$((($NJOB+1)*$PBS_NCPUS*$TASKS_THIS_ITERATION))

if [ $CUR_MAX -gt $CUR_TOT ]; then
    CUR_MAX=$CUR_TOT
fi

sed -n -e "${CUR_MIN},${CUR_MAX}p" $CMDS_PATH > ./JOB_PATH.txt
mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file ./JOB_PATH.txt --timeout 172800

$ECHO "All jobs finished, moving output..."

# Combine output
cd /scratch/ht96/nb9894/$TOTALJOBNAME/

cat ./*_mrr_done.csv >> $SAVEDIR/out_h2_mrr.csv
cat ./*_mkr_done.csv >> $SAVEDIR/out_h2_mkr.csv

# Delete loose files with
rm ./*_done.csv

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
    qsub -v NJOBS=$NJOBS,NJOB=$NJOB ./${JOBNAME}.sh
else
    $ECHO "Finished last job in sequence of $NJOBS jobs"
fi