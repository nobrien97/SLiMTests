#!/bin/bash -l
#PBS -P ht96
#PBS -l walltime=24:00:00
#PBS -l ncpus=1440
#PBS -l mem=2000GB
#PBS -l jobfs=1000GB
#PBS -l storage=scratch/ht96+gdata/ht96

ECHO=/bin/echo
JOBNAME=equalK_sweep/equalK_ODE/indPheno

mkdir -p /scratch/ht96/nb9894/$JOBNAME
mkdir -p /g/data/ht96/nb9894/$JOBNAME
mkdir -p $HOME/tests/$JOBNAME/done

cd $PBS_O_WORKDIR
SAVEDIR=/g/data/ht96/nb9894/$JOBNAME

$ECHO "Running nci-parallel..."
# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load nci-parallel/1.0.0a
export ncores_per_task=1
export ncores_per_numanode=12

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file ./cmds.txt --timeout 172800

$ECHO "All jobs finished, moving output..."

# Combine output into a single file
cd /scratch/ht96/nb9894/$JOBNAME

cat ./ode_* >> $SAVEDIR/ode_solutions_inds.csv

find -regex ".*[0-9]*_*[0-9].csv+" -delete
