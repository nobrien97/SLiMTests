#!/bin/bash -l
#PBS -P ht96
#PBS -l walltime=5:00:00
#PBS -l ncpus=5904
#PBS -l mem=10000GB
#PBS -l jobfs=6240GB
#PBS -l storage=scratch/ht96+gdata/ht96


# These variables are assumed to be set:
#   NJOBS is the total number of jobs in a sequence of jobs (defaults to 1)
#   NJOB is the number of the current job in the sequence (defaults to 0)

JOBNAME=getH2_mol_fix_lvls
JOBNAME_PREFIX=h2_mol_fix_lvls

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
    mkdir /scratch/ht96/nb9894/$JOBNAME_PREFIX/$JOBNAME
    mkdir /g/data/ht96/nb9894/$JOBNAME_PREFIX/$JOBNAME
    mkdir $HOME/tests/$JOBNAME_PREFIX/$JOBNAME/done

fi

# Quick terminate job sequence
if [ -f STOP_SEQUENCE ] ; then
    $ECHO  "Terminating sequence at job number $NJOB"
    exit 0
fi


cd $PBS_O_WORKDIR


SAVEDIR=/g/data/ht96/nb9894/$JOBNAME_PREFIX/$JOBNAME

# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load nci-parallel/1.0.0a
export ncores_per_task=1
export ncores_per_numanode=12

# Calculate the range of parameter combinations we are exploring this job

# CAUTION: may error if CUR_TOT is not a multiple of PBS_NCPUS - untested
CMDS_PATH=$HOME/tests/$JOBNAME_PREFIX/$JOBNAME/PBS/cmds.txt
CMD_LEN=$(cat $CMDS_PATH | wc -l)
CMD_MIN=$((($CMD_LEN/($NJOBS+1))*($NJOB) + 1))
CMD_MAX=$((($CMD_LEN/($NJOBS+1))*($NJOB+1)))
sed -n -e "${CMD_MIN},${CMD_MAX}p" $CMDS_PATH > ./JOB_PATH.txt


mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --dedicated --input-file ./JOB_PATH.txt --timeout 172800

# Combine output into a single file
cd /scratch/ht96/nb9894/$JOBNAME_PREFIX/$JOBNAME

cat ./*_done.csv >> $SAVEDIR/out_h2.csv

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
    qsub -v NJOBS=$NJOBS,NJOB=$NJOB ./$JOBNAME.sh
else
    $ECHO "Finished last job in sequence of $NJOBS jobs"
fi