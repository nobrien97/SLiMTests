#!/bin/bash -l
#PBS -P ht96
#PBS -l walltime=24:00:00
#PBS -l ncpus=4
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -l storage=scratch/ht96+gdata/ht96

# Job to filter data

cd $PBS_O_WORKDIR

module load R/4.0.0

Rscript ../R/make_figures.r