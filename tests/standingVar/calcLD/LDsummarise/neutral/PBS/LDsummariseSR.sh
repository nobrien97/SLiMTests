#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
FILENAME=${MODELINDEX}
JOBNAME=standingVar/calcLD/LDsummarise/neutral
TESTDIR=$HOME/tests/$JOBNAME

echo "Beginning run $MODELINDEX at $(date)"

RSCRIPTNAME1=$HOME/tests/standingVar/calcLD/LDsummarise/neutral/R/sumLD.R
RSCRIPTNAME2=$HOME/tests/standingVar/calcLD/LDsummarise/neutral/R/sumLDfreq.R


if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Run the model
echo "Calculating output for modelindex = $MODELINDEX...\n"

# Calculate stats for this model set
Rscript ${RSCRIPTNAME1} ${MODELINDEX}
Rscript ${RSCRIPTNAME2} ${MODELINDEX}

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX finished at $(date)!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
