#!/bin/bash -l
#PBS -P ht96
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l ncpus=728
#PBS -l mem=3500GB
#PBS -l jobfs=2800GB
#PBS -l storage=scratch/ht96+gdata/ht96
  
  
ECHO=/bin/echo
JOBNAME=newMotifs
FULLJOBNAME=$JOBNAME
#
# These variables are assumed to be set:
#   NJOBS is the total number of jobs in a sequence of jobs (defaults to 1)
#   NJOB is the number of the current job in the sequence (defaults to 0)
#   For this job, NJOBS should = 0
  
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
    mkdir -p /scratch/ht96/nb9894/$FULLJOBNAME
    mkdir -p /g/data/ht96/nb9894/$FULLJOBNAME
    mkdir -p $HOME/tests/$FULLJOBNAME/done
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
SAVEDIR=/g/data/ht96/nb9894/$FULLJOBNAME


# ========================================================================
# .... USER INSERTION OF EXECUTABLE LINE HERE 
# ========================================================================
# Make sure we're at the right place so we can find the bash script to run

$ECHO "Running nci-parallel..."
# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load nci-parallel/1.0.0a
export ncores_per_task=1
export ncores_per_numanode=13

# Calculate the range of parameter combinations we are exploring this job
# CAUTION: may error if CUR_TOT is not a multiple of PBS_NCPUS - untested
CMDS_PATH=$HOME/tests/$FULLJOBNAME/PBS/cmds.txt
CUR_TOT=$(cat $CMDS_PATH | wc -l)
CUR_MIN=$(($NJOB*$PBS_NCPUS+1))
CUR_MAX=$((($NJOB+1)*$PBS_NCPUS))

if [ $CUR_MAX -gt $CUR_TOT ]; then
    CUR_MAX=$CUR_TOT
fi

sed -n -e "${CUR_MIN},${CUR_MAX}p" $CMDS_PATH > ./JOB_PATH.txt

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file ./JOB_PATH.txt --timeout 172800

$ECHO "All jobs finished, moving output..."

# Combine output into a single file
cd /scratch/ht96/nb9894/$FULLJOBNAME/

cat ./slim_pos* >> $SAVEDIR/slim_pos.csv
cat ./slim_opt* >> $SAVEDIR/slim_opt.csv
cat ./slim_muts* >> $SAVEDIR/slim_muts.csv
cat ./slim_qg* >> $SAVEDIR/slim_qg.csv
cat ./slim_indPheno* >> $SAVEDIR/slim_indPheno.csv
cat ./slim_haplo_fix* >> $SAVEDIR/slim_haplo_fix.csv
# Remove haplo fix before we collect regular haplos - otherwise they get stuck to the end of the file
rm ./slim_haplo_fix*
cat ./slim_haplo* >> $SAVEDIR/slim_haplo.csv
cat ./slim_sampled_pheno* >> $SAVEDIR/slim_sampled_pheno.csv
cat ./slim_sampled_moltrait* >> $SAVEDIR/slim_sampled_moltrait.csv
cat ./slim_fx* >> $SAVEDIR/slim_fx.csv
cat ./slim_locusHo* >> $SAVEDIR/slim_locusHo.csv
cat ./slim_PMmat* >> $SAVEDIR/slim_PMmat.csv
cat ./slim_relPos* >> $SAVEDIR/slim_relPos.csv
cat ./slim_relVals* >> $SAVEDIR/slim_relVals.csv
cat ./slim_sharedmutfreqs* >> $SAVEDIR/slim_sharedmutfreqs.csv

# Save population state
mkdir -p $SAVEDIR/popstates
mv ./slim_popstate* $SAVEDIR/popstates

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
    qsub -v NJOBS=$NJOBS,NJOB=$NJOB ./$JOBNAME.sh
else
    $ECHO "Finished last job in sequence of $NJOBS jobs"
fi