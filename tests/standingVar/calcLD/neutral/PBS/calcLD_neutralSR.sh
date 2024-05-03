#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
SEED=$1
MODELINDEX=$2
FILENAME=${SEED}_${MODELINDEX}

BASENAME=standingVar/calcLD
JOBNAME=neutral
BASEDIR=$HOME/tests/$BASENAME
TESTDIR=$BASEDIR/$JOBNAME

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct seed from the file: put into array
SEED_FILE=$TESTDIR/R/standingVar_seeds.csv
SEED_NUM=($(awk "NR==$SEED" $SEED_FILE))

# Get correct model
MODEL_FILE=$TESTDIR/R/combos.csv
MODEL_NUM=($(awk "NR==$MODELINDEX" $MODEL_FILE))

# Run the model
echo "Running seed = $SEED...\n"
# If we have a K model, we need to disable molTraitFix by setting it to -1
if [[ "${MODEL_NUM[0]}" == "'K'" ]]
then
    $HOME/SLiM/slim -s ${SEED_NUM} -d modelindex=$MODELINDEX -d molTraitFix=-1 -d molTraitProps="c(0.25, 0.25, 0.25, 0.25)" -d modelType="'ODE'" -d rwide=${MODEL_NUM[1]} $TESTDIR/slim/calcLD_neutral.slim
else
    $HOME/SLiM/slim -s ${SEED_NUM} -d modelindex=$MODELINDEX -d molTraitFix="c(2,3)" -d molTraitProps="c(0.5, 0.5, 0.0, 0.0)" -d modelType="${MODEL_NUM[0]}" -d rwide=${MODEL_NUM[1]} $TESTDIR/slim/calcLD_neutral.slim
fi

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
