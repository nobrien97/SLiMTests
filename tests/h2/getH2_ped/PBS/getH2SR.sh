#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0

GDATAPATH=/g/data/ht96/nb9894/h2
SCRATCHPATH=/scratch/ht96/nb9894/h2/getH2_ped

# Two arguments: 
# run number (starts at 0)
# chunk (starts at 1): 
#        since we have heaps of files, we want to make sure that we don't write too many for /scratch/
#        so we chunk our data into sets so we can progressively combine files and delete them without
#        exceeding the limit.
RUN=$1
CHUNK=$2

if [ -f $HOME/tests/h2/getH2_ped/done/${RUN}_${CHUNK} ]; then
    echo "$RUN already done! Moving to next simulation."
    exit 0
fi

# Get the correct range of lines from the big files
## 203616 total models, 101 timepoints

## First get the min/max rows to read: pedigree is offset by 1 because header, haplotype doesn't have one
#MIN_HAPLO=$((1+$RUN*1000))
#MAX_HAPLO=$(((1+$RUN)*1000))
MIN_PED=$((1+$RUN*500))
MAX_PED=$(((1+$RUN)*500))

# Save files
#tail -n "+$MIN_HAPLO" $GDATAPATH/slim_haplo.csv | head -n "$(($MAX_HAPLO-$MIN_HAPLO+1))" > slim_haplo_sbst_$RUN.csv
tail -n "+$MIN_PED" $GDATAPATH/slim_pedigree.csv | head -n "$(($MAX_PED-$MIN_PED+1))" > slim_ped_sbst_$RUN.csv
tail -n "+$(($RUN+1))" $GDATAPATH/slim_sampled_pheno.csv | head -n 1 > slim_pheno_sbst_$RUN.csv

RSCRIPTNAME=$HOME/tests/h2/getH2_ped/R/calcH2.R

echo "Running RUNindex = $RUN...\n"
Rscript ${RSCRIPTNAME} ${RUN} ${CHUNK}


# Create file to show what we've already done if we get interrupted
touch $HOME/tests/h2/getH2_ped/done/${RUN}_${CHUNK}

# Check if we're the last in a chunk, if we are we need to do some cleanup, otherwise we can continue
# Chunks should consist of 2016 files
if [ $(ls $HOME/tests/h2/getH2_ped/done/*_${CHUNK} | wc -l) == 2016 ]; then
    echo "Chunk $CHUNK done, combining chunk files and cleaning up..."
    cat $SCRATCHPATH/*_${CHUNK}.csv >> $SCRATCHPATH/out_h2_${CHUNK}_done.csv
    rm $SCRATCHPATH/*_${CHUNK}.csv
fi


DURATION=$SECONDS
echo "Run RUNindex = $RUN, finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."
