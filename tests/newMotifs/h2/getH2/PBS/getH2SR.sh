#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0
SUBJOBNAME=newMotifs
JOBNAME=h2/getH2
TOTALJOBNAME=$SUBJOBNAME/$JOBNAME



# Two arguments: 
# run number (starts at 0)
# chunk (starts at 1): 
#        since we have heaps of files, we want to make sure that we don't write too many for /scratch/
#        so we chunk our data into sets so we can progressively combine files and delete them without
#        exceeding the limit.
RUN=$1
CHUNK=$2

if [ -f $HOME/tests/${TOTALJOBNAME}/done/${RUN}_* ]; then
    echo "$RUN already done! Moving to next simulation."
    exit 0
fi

SCRATCHPATH=/scratch/ht96/nb9894/${TOTALJOBNAME}
DATAPATH=/g/data/ht96/nb9894/${SUBJOBNAME}
#DATAPATH=./

fix_haplo=$(tail -n "+${RUN}" $DATAPATH/slim_haplo_fix.csv | head -n 1)

# Save subset files to work on
echo $fix_haplo > slim_haplo_fix_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_haplo.csv | head -n 1 > slim_haplo_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_relVals.csv | head -n 1 > slim_relVals_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_relPos.csv | head -n 1 > slim_relPos_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_sampled_moltrait.csv | head -n 1 > slim_moltrait_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_sampled_pheno.csv | head -n 1 > slim_pheno_sbst_$RUN.csv


RSCRIPTNAME=$HOME/tests/${TOTALJOBNAME}/R/calcH2.R

echo "Running RUNindex = $RUN...\n"
Rscript ${RSCRIPTNAME} ${RUN} ${CHUNK}

# Create file to show what we've already done if we get interrupted
touch $HOME/tests/${TOTALJOBNAME}/done/${RUN}_${CHUNK}

# Check if we're the last in a chunk, if we are we need to do some cleanup, otherwise we can continue
# Chunks should consist of 1440 files
# if [ $(ls $HOME/tests/${TOTALJOBNAME}/done/*_${CHUNK} | wc -l) -ge 1440 ]; then
#     echo "Chunk $CHUNK done, combining chunk files and cleaning up..."
#     cat $SCRATCHPATH/*_${CHUNK}_mrr.csv >> $SCRATCHPATH/out_h2_${CHUNK}_mrr_done.csv
#     cat $SCRATCHPATH/*_${CHUNK}_mkr.csv >> $SCRATCHPATH/out_h2_${CHUNK}_mkr_done.csv
#     rm $SCRATCHPATH/*_${CHUNK}_mrr.csv
#     rm $SCRATCHPATH/*_${CHUNK}_mkr.csv
# fi

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished at $(date)!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."
