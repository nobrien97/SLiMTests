#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
FILENAME=${MODELINDEX}
JOBNAME=calcMutationStats
TESTDIR=$HOME/tests/$JOBNAME

RSCRIPTNAME=$HOME/tests/standingVar/analysis/R/calcMutationStats.R

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct modelindex from the file: put into array
MODEL_FILE=$TESTDIR/R/combos.csv
MODEL_NUM=($(awk "NR==$MODELINDEX" $MODEL_FILE))

# Run the model
echo "Calculating output for modelindex = $MODELINDEX...\n"

# Calculate stats for this model set
Rscript ${RSCRIPTNAME} ${MODELINDEX}

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
