#!/bin/bash -l
#PBS -P ht96
#PBS -q normalsr
#PBS -l walltime=24:00:00
#PBS -l ncpus=104
#PBS -l mem=500GB
#PBS -l jobfs=400GB
#PBS -l storage=scratch/ht96+gdata/ht96

JOBNAME=newMotifs/paper1/ruggedness/sim_range/pca_landscape
TESTDIR=$HOME/tests/$JOBNAME

# Calculate landscape
module load R/4.0.0
RSCRIPTNAME=$TESTDIR/R/pca_landscape.R
Rscript ${RSCRIPTNAME}

