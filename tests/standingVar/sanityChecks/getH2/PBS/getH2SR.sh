#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0
SUBJOBNAME=standingVar/sanityChecks
JOBNAME=getH2
TOTALJOBNAME=$SUBJOBNAME/$JOBNAME



# one arguments: 
# run number (starts at 0)
RUN=$1

if [ -f $HOME/tests/${TOTALJOBNAME}/done/${RUN}_* ]; then
    echo "$RUN already done! Moving to next simulation."
    exit 0
fi

SCRATCHPATH=/scratch/ht96/nb9894/${SUBJOBNAME}
DATAPATH=/g/data/ht96/nb9894/${SUBJOBNAME}

# Save subset files to work on
tail -n "+${RUN}" $DATAPATH/slim_haplo_fix.csv | head -n 1 > slim_haplo_fix_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_haplo.csv | head -n 1 > slim_haplo_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_relVals.csv | head -n 1 > slim_relVals_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_relPos.csv | head -n 1 > slim_relPos_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_sampled_moltrait.csv | head -n 1 > slim_moltrait_sbst_$RUN.csv
tail -n "+${RUN}" $DATAPATH/slim_sampled_pheno.csv | head -n 1 > slim_pheno_sbst_$RUN.csv

RSCRIPTNAME=$HOME/tests/${TOTALJOBNAME}/R/calcH2.R

echo "Running RUNindex = $RUN...\n"
Rscript ${RSCRIPTNAME} ${RUN}

# Create file to show what we've already done if we get interrupted
touch $HOME/tests/${TOTALJOBNAME}/done/${RUN}

DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seed = $SEED finished at $(date)!"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."
