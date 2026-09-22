#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODEL=$1
FILENAME=${MODEL}
JOBNAME=newMotifs/paper1/ruggedness/umap
TESTDIR=$HOME/tests/$JOBNAME

echo "Beginning run model = $MODEL at $(date)"

RSCRIPTNAME=$TESTDIR/R/ruggedness_umap.R

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Run the model
echo "Calculating output for modelindex = $MODELINDEX...\n"

# Calculate stats for this model set
Rscript ${RSCRIPTNAME} ${MODEL}

DURATION=$SECONDS
echo "Run modelindex = $MODEL finished at $(date)!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
