#!/bin/bash -l

module load R/4.2.2

cd $PBS_JOBFS
SECONDS=0

RUN=$1
FILENAME=${RUN}

JOBNAME=equalK_sweep/equalK_ODE/indPheno
TESTDIR=$HOME/tests/$JOBNAME

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi


SCRATCHPATH=/scratch/ht96/nb9894/$JOBNAME

RSCRIPTNAME=${TESTDIR}/R/ODESolution.R

echo "Running RUNindex = $RUN...\n"
Rscript ${RSCRIPTNAME} ${RUN}

# Move output
mv ode_${RUN}.csv ${SCRATCHPATH}

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${RUN}

DURATION=$SECONDS
echo "Run RUNindex = $RUN, finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."
