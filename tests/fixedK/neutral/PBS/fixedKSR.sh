#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
SEED=$2
FILENAME=${MODELINDEX}_${SEED}
JOBNAME=fixedK/neutral
TESTDIR=$HOME/tests/$JOBNAME

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct modelindex from the file: put into array
MODEL_FILE=$TESTDIR/R/combos.csv
MODEL_NUM=($(awk "NR==$MODELINDEX" $MODEL_FILE))
SEED_FILE=$TESTDIR/R/${JOBNAME}_seeds.csv
SEED_NUM=($(awk "NR==$SEED" $SEED_FILE))

MODEL=indTrack_net.slim

# if [[ "${MODEL_NUM[0]}" != "Additive" ]]; then
#    MODEL=indTrack_net.slim
#    MODEL_NUM[1]=$((MODEL_NUM[1]-2))
# fi

echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"
# Run the model - set width to 0 for neutral simulation 
$HOME/SLiM/slim -s ${SEED_NUM} -d modelindex=$MODELINDEX -d nloci=${MODEL_NUM[0]} -d width=0 $TESTDIR/slim/$MODEL

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
