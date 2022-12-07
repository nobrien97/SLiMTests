#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0

JOBNAME=getH2_mol_fix_lvls
JOBNAME_PREFIX=h2_mol_fix_lvls


GDATAPATH=/g/data/ht96/nb9894/$JOBNAME_PREFIX
SCRATCHPATH=/scratch/ht96/nb9894/$JOBNAME_PREFIX/$JOBNAME

# Three arguments: 
# run number (starts at 0)
# chunk (starts at 1): 
#        since we have heaps of files, we want to make sure that we don't write too many for /scratch/
#        so we chunk our data into sets so we can progressively combine files and delete them without
#        exceeding the limit.
RUN=$1
CHUNK=$2

if [ -f $HOME/tests/$JOBNAME_PREFIX/$JOBNAME/done/${RUN}_* ]; then
    echo "$RUN already done! Moving to next simulation."
    exit 0
fi

# Get the correct range of lines from the big files
## 5904 total models

# Choose the correct file

N_SAMP=2000

## First get the min/max rows to read: offset by HAPLO_NUM
MIN_HAPLO=$((1+$RUN*$N_SAMP))
MAX_HAPLO=$(((1+$RUN)*$N_SAMP))


# Save files
tail -n "+$MIN_HAPLO" $GDATAPATH/slim_haplo.csv | head -n "$(($MAX_HAPLO-$MIN_HAPLO+1))" > slim_haplo_sbst_$RUN.csv
# tail -n "+$MIN_PED" $GDATAPATH/slim_pedigree.csv | head -n "$(($MAX_PED-$MIN_PED+1))" > slim_ped_sbst_$RUN.csv
tail -n "+$(($RUN+1))" $GDATAPATH/slim_sampled_pheno.csv | head -n 1 > slim_pheno_sbst_$RUN.csv

RSCRIPTNAME=$HOME/tests/$JOBNAME_PREFIX/$JOBNAME/R/calcH2.R

echo "Running RUNindex = $RUN...\n"
Rscript ${RSCRIPTNAME} ${RUN} ${CHUNK}

# Create file to show what we've already done if we get interrupted
touch $HOME/tests/$JOBNAME_PREFIX/$JOBNAME/done/${RUN}_${CHUNK}

# Check if we're the last in a chunk, if we are we need to do some cleanup, otherwise we can continue
# Chunks should consist of 1968 files
if [ $(ls $HOME/tests/$JOBNAME_PREFIX/$JOBNAME/done/*_${CHUNK} | wc -l) == 1968 ]; then
    echo "Chunk $CHUNK done, combining chunk files and cleaning up..."
    cat $SCRATCHPATH/*_${CHUNK}.csv >> $SCRATCHPATH/out_h2_${CHUNK}_done.csv
    rm $SCRATCHPATH/*_${CHUNK}.csv
fi


DURATION=$SECONDS
echo "Run RUNindex = $RUN, finished!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."
