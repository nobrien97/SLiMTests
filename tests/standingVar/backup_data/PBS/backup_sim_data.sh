#!/bin/bash -l
#PBS -P ht96
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -l ncpus=1
#PBS -l mem=190GB
#PBS -l storage=scratch/ht96+gdata/ht96

# copy across data
SIMDATA_FILEPATH=/g/data/ht96/nb9894/standingVar
scp -i ~/.ssh/rdmgadi_id_rsa $SIMDATA_FILEPATH/standingVarMuts.db $SIMDATA_FILEPATH/*.csv uqnobri4@data.qriscloud.org.au:/QRISdata/Q4117/standingVar

