#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
SEED=$2
FILENAME=${MODELINDEX}_${SEED}
BASENAME=newMotifs
JOBNAME=mutVar
BASEDIR=$HOME/tests/$BASENAME
TESTDIR=$BASEDIR/$JOBNAME

COMBOSDIR=$HOME/tests/$BASENAME
ORIGJOBPATH=$HOME/tests/newMotifs/randomisedStarts
DATAPATH=/g/data/ht96/nb9894/newMotifs/randomisedStarts

echo "Beginning run modelindex = $MODELINDEX, seed = $SEED at $(date)"

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct modelindex from the file: put into array
MODEL_FILE=$COMBOSDIR/R/combos.csv
MODEL_NUM=($(awk "NR==$MODELINDEX" $MODEL_FILE))
SEED_FILE=$ORIGJOBPATH/R/${BASENAME}_seeds.csv
SEED_NUM=($(awk "NR==$SEED" $SEED_FILE))

# Extract data from /g/data
# Ensure that all three only have one row, if not choose the first (can happen if a simulation is repeated)
awk -F',' -v a="$MODELINDEX" -v b="$SEED_NUM" '{if ($1 == a && $2 == b){print;exit;}}' $DATAPATH/slim_pos.csv > slim_pos${MODELINDEX}_${SEED}.csv
awk -F',' -v a="$MODELINDEX" -v b="$SEED_NUM" '{if ($1 == b && $2 == a){print;exit;}}' $DATAPATH/slim_opt.csv > slim_opt${MODELINDEX}_${SEED}.csv

# Run the model
echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"

$HOME/SLiM/chp3/slim -s ${SEED_NUM} -d modelindex=$MODELINDEX -d inputSeed=$SEED -d inputModel=$MODELINDEX -d modelType="${MODEL_NUM[0]}" $TESTDIR/slim/calcMutVar.slim

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished at $(date)!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
