#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
SEED=$1
MODELINDEX=$2
FILENAME=${SEED}_${MODELINDEX}

if [ -f $HOME/tests/h2_ext/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

RSCRIPTNAME=$HOME/tests/h2_ext/R/h2_singleshot.R

echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"
Rscript ${RSCRIPTNAME} ${SEED} ${MODELINDEX}


DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $HOME/tests/h2_ext/done/${FILENAME}
