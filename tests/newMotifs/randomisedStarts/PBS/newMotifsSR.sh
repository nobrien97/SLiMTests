#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
SEED=$2
FILENAME=${MODELINDEX}_${SEED}
JOBNAME=newMotifs
FULLJOBNAME=$JOBNAME/randomisedStarts
TESTDIR=$HOME/tests/$FULLJOBNAME
COMBOSDIR=$HOME/tests/$JOBNAME

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct modelindex from the file: put into array
MODEL_FILE=$COMBOSDIR/R/combos.csv
MODEL_NUM=($(awk "NR==$MODELINDEX" $MODEL_FILE))
SEED_FILE=$TESTDIR/R/${JOBNAME}_seeds.csv
SEED_NUM=($(awk "NR==$SEED" $SEED_FILE))

# Run the model
echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"
$HOME/SLiM/slim -s ${SEED_NUM} -d modelindex=$MODELINDEX -d molTraitFix=-1 -d modelType=${MODEL_NUM[0]} -d rwide=${MODEL_NUM[1]} $TESTDIR/slim/baseScript.slim

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
