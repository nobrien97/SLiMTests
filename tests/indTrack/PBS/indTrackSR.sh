#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
FILENAME=${MODELINDEX}
JOBNAME=indTrack
TESTDIR=$HOME/tests/$JOBNAME

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct modelindex from the file: put into array
MODEL_FILE=$TESTDIR/R/combos.csv
MODEL_NUM=($(awk "NR==$MODELINDEX" $MODEL_FILE))

MODEL=indTrack_add.slim
if [[ "${MODEL_NUM[0]}" != "Additive" ]]; then
    MODEL=indTrack_net.slim
    MODEL_NUM[1]=$((MODEL_NUM[1]-2))
fi

echo "Running modelindex = $MODELINDEX...\n"
# Run the model
$HOME/SLiM/slim -s ${MODEL_NUM[3]} -d modelindex=$MODELINDEX -d nloci=${MODEL_NUM[1]} -d locisigma=${MODEL_NUM[2]} $TESTDIR/slim/$MODEL


DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
