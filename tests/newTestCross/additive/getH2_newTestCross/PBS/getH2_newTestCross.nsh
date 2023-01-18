#!/sw7/RCC/NimrodG/embedded-1.11.0/bin/nimexec

#Fix up Account String
#PBS -A UQ-SCI-BiolSci
#6 nodes, 24 cores, 120GB per node
#PBS -l select=6:ncpus=24:mem=120GB:ompthreads=1 
#PBS -l walltime=48:00:00
#PBS -N SLIM_NIM


# =============================================================================
# Tell Nimrod what range of parameter values you want to use 
# =============================================================================
# There are additional directives for Nimrod to interpret with #NIM at the start of each line.
# Tell Nimrod to use this as the shell for the job proper when it has parsed this file.
#NIM shebang /bin/bash

#The parameter for the row in the input file to get
#NIM parameter RUN integer range from 0 to 11807 step 1

#The chunk - only splitting this one in four since it's relatively small
#NIM parameter SPLIT integer range from 0 to 3 step 1

# =============================================================================
# Code for individual run
# =============================================================================
cd $PBS_O_WORKDIR

# Just checking that something did not go wrong with assignment of the J values.
if [ -z "${NIMROD_VAR_RUN}" ]; then
        echo "\$NIMROD_VAR_RUN isn't set, cannot continue..." 
        exit 2
fi

if [ -z "${NIMROD_VAR_SPLIT}" ]; then
        echo "\$NIMROD_VAR_SPLIT isn't set, cannot continue..." 
        exit 2
fi

#Where you submit this job from will be the value of $PBS_O_WORKDIR
echo "PBS_O_WORKDIR is ${PBS_O_WORKDIR}" 


module load R/4.0.0

cd $TMPDIR
SECONDS=0

JOBNAME=getH2_newTestCross
JOBNAME_PREFIX=newTestCross/additive


SCRATCHPATH=/scratch/user/uqnobri4/$JOBNAME_PREFIX

# Three arguments: 
# run number (starts at 0)
# haplo_num: since the haplotypes file is huge, I've split it into some smaller files that will be faster
#        to process. This means we need to keep track of what file to open, and also use that info to offset
#        our line selection, which we use HAPLO_NUM for.
RUN=$NIMROD_VAR_RUN
HAPLO_NUM=$NIMROD_VAR_SPLIT

if [ -f $HOME/tests/$JOBNAME_PREFIX/$JOBNAME/done/${RUN} ]; then
    echo "$RUN already done! Moving to next simulation."
    exit 0
fi

# Get the correct range of lines from the big files
## 11808 total models 

# Choose the correct file
HAPLO_FILE=x0$HAPLO_NUM

N_SAMP=2000

SUBSET_LEN=$((11808/4)) # Kinda hacky, can I directly get CMD_LEN from the main script?

## First get the min/max rows to read: offset by HAPLO_NUM
MIN_HAPLO=$(((1+$RUN*$N_SAMP)-$HAPLO_NUM*$SUBSET_LEN*$N_SAMP))
MAX_HAPLO=$((((1+$RUN)*$N_SAMP)-$HAPLO_NUM*$SUBSET_LEN*$N_SAMP))

# Save files
tail -n "+$MIN_HAPLO" $SCRATCHPATH/$HAPLO_FILE | head -n "$(($MAX_HAPLO-$MIN_HAPLO+1))" > slim_haplo_sbst_$RUN.csv
tail -n "+$(($RUN+1))" $SCRATCHPATH/slim_sampled_pheno.csv | head -n 1 > slim_pheno_sbst_$RUN.csv

RSCRIPTNAME=$HOME/tests/$JOBNAME_PREFIX/$JOBNAME/R/calcH2.R

echo "Running RUNindex = $RUN...\n"
Rscript ${RSCRIPTNAME} ${RUN}

# Create file to show what we've already done if we get interrupted
touch $HOME/tests/$JOBNAME_PREFIX/$JOBNAME/done/${RUN}


DURATION=$SECONDS
echo "Run RUNindex = $RUN, finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."
