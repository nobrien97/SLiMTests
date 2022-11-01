#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0

GDATAPATH=/g/data/ht96/nb9894/h2

# Only one argument in this case: the run number (starts at 0) 
RUN=$1

if [ -f $HOME/tests/h2/getH2/done/${RUN} ]; then
    echo "$RUN already done! Moving to next simulation."
    exit 0
fi

# Get the correct range of lines from the big files
## 203616 total models, 101 timepoints

## First get the min/max rows to read: pedigree is offset by 1 because header, haplotype doesn't have one
MIN_HAPLO=$((1+$RUN*1000))
MAX_HAPLO=$(((1+$RUN)*1000))
MIN_PED=$(((1+$RUN*500)+1))
MAX_PED=$((((1+$RUN)*500)+1))

# Save files
tail -n "+$MIN_HAPLO" $GDATAPATH/slim_haplo.csv | head -n "$(($MAX_HAPLO-$MIN_HAPLO+1))" > slim_haplo_sbst_$RUN.csv
tail -n "+$MIN_PED" $GDATAPATH/slim_pedigree.csv | head -n "$(($MAX_PED-$MIN_PED+1))" > slim_ped_sbst_$RUN.csv

RSCRIPTNAME=$HOME/tests/h2/getH2/R/h2_singleshot.R

echo "Running RUNindex = $RUN...\n"
Rscript ${RSCRIPTNAME} ${RUN}


DURATION=$SECONDS
echo "Run RUNindex = $RUN, finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."

# Create file to show what we've already done if we get interrupted
touch $HOME/tests/h2/getH2/done/${RUN}
