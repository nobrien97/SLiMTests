#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
SEED=$2
POPFILE=$3
FILENAME=${MODELINDEX}_${SEED}
JOBNAME=fixedK
TESTDIR=$HOME/tests/$JOBNAME/calcHet

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi


MODEL=calcHetFromPopState_add.slim
NLOCI=2
# If we are asking for more than 2 loci, we're running a network model and need to switch files
if [[ ${MODELINDEX} == 2 ]]; then
   MODEL=calcHetFromPopState_net.slim
   NLOCI=4
fi

echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"
# Run the model
$HOME/SLiM/slim -s ${SEED} -d modelindex=$MODELINDEX -d nloci=${NLOCI} -d popPath=\"${POPFILE}\" $TESTDIR/slim/$MODEL

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
