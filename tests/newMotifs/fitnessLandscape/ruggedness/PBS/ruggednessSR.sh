#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
FILENAME=${MODELINDEX}
JOBNAME=newMotifs/fitnessLandscape/ruggedness
TESTDIR=$HOME/tests/$JOBNAME

echo "Beginning run modelindex = $MODELINDEX at $(date)"

RSCRIPTNAME=$HOME/tests/newMotifs/fitnessLandscape/ruggedness/R/nosil_permolcomp_parallel.R

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Calculate stats for this model set
Rscript ${RSCRIPTNAME} ${MODELINDEX}

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX finished in $(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds!"

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
