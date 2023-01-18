#!/sw7/RCC/NimrodG/embedded-1.11.0/bin/nimexec

#Fix up Account String
#PBS -A UQ-SCI-BiolSci
#6 nodes, 24 cores, 120GB per node
#PBS -l select=6:ncpus=24:mem=120GB:ompthreads=1 
#PBS -l walltime=48:00:00
#PBS -N SLIM_NIM

# There are additional directives for Nimrod to interpret with #NIM at the start of each line.
# Tell Nimrod to use this as the shell for the job proper when it has parsed this file.
#NIM shebang /bin/bash

COMBO_ROWS=6
SEED_ROWS=48

# =============================================================================
# Tell Nimrod what range of parameter values you want to use 
# =============================================================================

#The parameters for the latin squares are rows in the input file.
#NIM parameter LS integer range from 1 to ${COMBO_ROWS} step 1

#Repeat 48 times for each hypercube with a different SEED value
#NIM parameter SEED integer range from 1 to ${SEED_ROWS} step 1


# Just checking that something did not go wrong with assignment of the J values.
if [ -z "${NIMROD_VAR_LS}" ]; then
        echo "\$NIMROD_VAR_LS isn't set, cannot continue..." 
        exit 2
fi

if [ -z "${NIMROD_VAR_SEED}" ]; then
        echo "\$NIMROD_VAR_SEED isn't set, cannot continue..." 
        exit 2
fi

#Where you submit this job from will be the value of $PBS_O_WORKDIR
echo "PBS_O_WORKDIR is ${PBS_O_WORKDIR}" 
#Everything you need should be located relative to PBS_O_WORKDIR, or else a full path

#Set the cd to TMPDIR for writing SLiM output
cd ${TMPDIR}

# Always run the entire parameter range cause nimrod can do them in any order.

SECONDS=0

# Load the right gcc so we have access to libstdc++6

module load gcc/8.5.0

# Rename the first and second arguments passed to this single shot script for clarity 
SEED=$NIMROD_VAR_SEED
MODELINDEX=$NIMROD_VAR_LS
FILENAME=${SEED}_${MODELINDEX}
JOBNAME=newTestCross
TESTDIR=$HOME/tests/$JOBNAME/additive

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct seed from the file: no header on this one, no need to adjust
SEED_FILE=$TESTDIR/R/${JOBNAME}_seeds.csv
SEED_NUM=$(awk "NR==$SEED" $SEED_FILE)

# Get the correct modelindex from the file: put into array
MODEL_FILE=$TESTDIR/R/combos.csv
MODEL_NUM=($(awk "NR==$MODELINDEX" $MODEL_FILE))

MODEL=hsfs_additive.slim

echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"
# Run the model
$HOME/slim -s $SEED_NUM -d modelindex=$MODELINDEX -d nloci=${MODEL_NUM[0]} -d locisigma=${MODEL_NUM[1]} $TESTDIR/slim/$MODEL


DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}

