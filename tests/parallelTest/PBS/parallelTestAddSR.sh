#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
SEED=$2
FILENAME=${MODELINDEX}_${SEED}
JOBNAME=parallelTest
TESTDIR=$HOME/tests/$JOBNAME

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct modelindex from the file: put into array
SEED_FILE=$TESTDIR/R/${JOBNAME}_seeds.csv
SEED_NUM=($(awk "NR==$SEED" $SEED_FILE))

echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"
# Run the model
$HOME/SLiM/slim_multi -s ${SEED_NUM} -d modelType='"Add"' -d modelindex=$MODELINDEX -maxThreads ${ncores_per_task} $TESTDIR/slim/baseScript.slim

DURATION=$SECONDS
echo "Run maxThreads = $ncores_per_task, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
