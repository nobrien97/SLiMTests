#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
SEED=$1
MODELINDEX=$2

FILENAME=${SEED}_${MODELINDEX}
JOBNAME=h2_mol_fix_lvls
TESTDIR=$HOME/tests/$JOBNAME

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct seed from the file: no header on this one, no need to adjust
SEED_FILE=$TESTDIR/R/${JOBNAME}_seeds.csv
SEED_NUM=$(awk "NR==$SEED" $SEED_FILE)

MODEL=hsfs_network.slim
# Get the correct molecular trait value based on modelindex (0.5, 1, 1.5, or 2)
MOLVAL=$(echo "(${MODELINDEX}+1)*0.5" | bc -l)

echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"
# Run the model
$HOME/SLiM/slim -s $SEED_NUM -d modelindex=$MODELINDEX -d molTraitFix=${MODELINDEX} -d molVal=${MOLVAL} $TESTDIR/slim/$MODEL


DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}