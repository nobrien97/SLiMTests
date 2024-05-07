#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
SEED=$2
FILENAME=${MODELINDEX}_${SEED}
BASENAME=standingVar
JOBNAME=calcLD
BASEDIR=$HOME/tests/$BASENAME
TESTDIR=$BASEDIR/$JOBNAME
DATAPATH=/g/data/ht96/nb9894/standingVar

echo "Beginning run modelindex = $MODELINDEX, seed = $SEED at $(date)"

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Get the correct modelindex from the file: put into array
MODEL_FILE=$BASEDIR/R/combos.csv
MODEL_NUM=($(awk "NR==$MODELINDEX" $MODEL_FILE))
SEED_FILE=$BASEDIR/R/${BASENAME}_seeds.csv
SEED_NUM=($(awk "NR==$SEED" $SEED_FILE))

# Extract data from /g/data
# Ensure that all three only have one row, if not choose the first (can happen if a simulation is repeated)
awk -F',' -v a="$MODELINDEX" -v b="$SEED_NUM" '{if ($1 == a && $2 == b){print;exit;}}' $DATAPATH/slim_pos.csv > slim_pos${SEED_NUM}_${MODELINDEX}.csv
awk -F',' -v a="$MODELINDEX" -v b="$SEED_NUM" '{if ($1 == a && $2 == b){print;exit;}}' $DATAPATH/slim_dict.csv > slim_dict${SEED_NUM}_${MODELINDEX}.csv
awk -F',' -v a="$MODELINDEX" -v b="$SEED_NUM" '{if ($1 == b && $2 == a){print;exit;}}' $DATAPATH/slim_opt.csv > slim_opt${SEED_NUM}_${MODELINDEX}.csv

# Run the model
echo "Running modelindex = $MODELINDEX, seed = $SEED...\n"


# If we have a K model, we need to disable molTraitFix by setting it to -1
if [[ "${MODEL_NUM[3]}" == "'K'" ]]
then
    $HOME/SLiM/slim -s ${SEED_NUM} -d modelindex=$MODELINDEX -d inputSeed=${SEED_NUM} -d inputModel=$MODELINDEX -d molTraitFix=-1 -d molTraitProps="c(0.25, 0.25, 0.25, 0.25)" -d nloci=${MODEL_NUM[0]} -d locisigma=${MODEL_NUM[1]} -d rwide=${MODEL_NUM[2]} -d modelType="'ODE'" $TESTDIR/slim/calcLD.slim
else
    $HOME/SLiM/slim -s ${SEED_NUM} -d modelindex=$MODELINDEX -d inputSeed=${SEED_NUM} -d inputModel=$MODELINDEX -d molTraitFix="c(2,3)" -d molTraitProps="c(0.5, 0.5, 0.0, 0.0)" -d nloci=${MODEL_NUM[0]} -d locisigma=${MODEL_NUM[1]} -d rwide=${MODEL_NUM[2]} -d modelType="${MODEL_NUM[3]}" $TESTDIR/slim/calcLD.slim
fi

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished at $(date)!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}
