#!/bin/bash -l

module load python3/3.12.1

export PYTHONPATH=/g/data/ht96/nb9894/py_libs/lib/python3.12/site-packages:$PYTHONPATH

# Activate virtual environment
source /g/data/ht96/nb9894/py_libs/virtual/environment/ruggedness_3.12.1/bin/activate

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODEL=$1
FILENAME=${MODEL}
JOBNAME=newMotifs/paper1/ruggedness/rhvae
TESTDIR=$HOME/tests/$JOBNAME

echo "Beginning run model = $MODEL at $(date)"

PYSCRIPTNAME=$TESTDIR/py/ruggedness_rhvae.py

if [ -f $TESTDIR/done/${FILENAME} ]; then
    echo "$FILENAME already done! Moving to next simulation."
    exit 0
fi

# Run the model
echo "Calculating output for modelindex = $MODELINDEX...\n"

# Calculate stats for this model set
python3 ${PYSCRIPTNAME} ${MODEL}

DURATION=$SECONDS
echo "Run modelindex = $MODEL finished at $(date)!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $TESTDIR/done/${FILENAME}

# Deactivate venv
deactivate