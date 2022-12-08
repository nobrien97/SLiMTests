#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
SEED=$1
MODELINDEX=$2
FILENAME=${SEED}_${MODELINDEX}
JOBNAME=newTestCross
TESTDIR=$HOME/tests/$JOBNAME

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

MODEL=hsfs_network.slim

echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"
# Run the model
$HOME/SLiM/slim -s $SEED_NUM -d modelindex=$MODELINDEX -d molTraitFix=${MODEL_NUM[0]} -d nloci=${MODEL_NUM[1]} -d locisigma=${MODEL_NUM[2]} $TESTDIR/slim/$MODEL


DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
