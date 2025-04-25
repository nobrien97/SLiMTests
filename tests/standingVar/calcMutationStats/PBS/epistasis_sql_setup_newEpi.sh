#!/bin/bash -l
#PBS -P ht96
#PBS -q hugemem
#PBS -l walltime=48:00:00
#PBS -l ncpus=4
#PBS -l mem=1470GB
#PBS -l iointensive=1
#PBS -l storage=scratch/ht96+gdata/ht96

# Purpose: Sets up a sql database for mutation data

# Path to updated version of sqlite
SQLITE3=~/Tools/sqlite/sqlite3
SAVEDIR=/g/data/ht96/nb9894/standingVar/calcMutationStats
RUNDIR=/iointensive

cd $PBS_O_WORKDIR

# Loads data into data frame and attaches indices
cat ./load_data_newEpi.sql | $SQLITE3 $RUNDIR/epistasis.db 2>/dev/null
cat ./add_indices_newEpi.sql | $SQLITE3 $RUNDIR/epistasis.db 2>/dev/null

mv $RUNDIR/epistasis.db* $SAVEDIR