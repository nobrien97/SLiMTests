#!/bin/bash -l
#SBATCH --account a_ortiz_barrientos_coe
#SBATCH --ntasks=14400
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --constraint=epyc3
#SBATCH --batch=epyc3
#SBATCH --partition=general

JOBNAME=standingVar/analysis
ECHO=/bin/echo

$ECHO "Creating outputs folders..."

# Make output folder
mkdir /scratch/user/uqnobri4/${JOBNAME}
mkdir /QRISdata/Q4117/${JOBNAME}
mkdir $HOME/tests/${JOBNAME}/done

SAVEDIR=/QRISdata/Q4117/${JOBNAME}

# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load r/4.2.1

export ncores_per_task=1
export ncores_per_numanode=12

# Calculate the range of parameter combinations we are exploring this job

# CAUTION: may error if CUR_TOT is not a multiple of PBS_NCPUS - untested
CMDS_PATH=$HOME/tests/${JOBNAME}/PBS/cmds.txt
CMD_LEN=$(cat $CMDS_PATH | wc -l)
CMD_MIN=$((($CMD_LEN/($NJOBS+1))*($NJOB) + 1))
CMD_MAX=$((($CMD_LEN/($NJOBS+1))*($NJOB+1)))
sed -n -e "${CMD_MIN},${CMD_MAX}p" $CMDS_PATH > ./JOB_PATH.txt


srun --ntasks=14400 --export=ALL 
mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --dedicated --input-file ./JOB_PATH.txt --timeout 172800

# Combine output into a single file
cd /scratch/user/uqnobri4/${JOBNAME}

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
    qsub -v NJOBS=$NJOBS,NJOB=$NJOB ./getH2_hsfs_nloci_split.sh
else
    $ECHO "Finished last job in sequence of $NJOBS jobs"
fi