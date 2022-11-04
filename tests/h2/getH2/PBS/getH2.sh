#!/bin/bash -l
#PBS -P ht96
#PBS -l walltime=5:00:00
#PBS -l ncpus=10368
#PBS -l mem=41040GB
#PBS -l jobfs=86400GB
#PBS -l storage=scratch/ht96+gdata/ht96
  
# create our folders
cd $PBS_O_WORKDIR

# Make output folder
mkdir /scratch/ht96/nb9894/h2/getH2
mkdir /g/data/ht96/nb9894/h2/getH2
mkdir $HOME/tests/h2/getH2/done

SAVEDIR=/g/data/ht96/nb9894/h2/getH2

# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load nci-parallel/1.0.0a
export ncores_per_task=1
export ncores_per_numanode=12

# Calculate the range of parameter combinations we are exploring this job
# CAUTION: may error if CUR_TOT is not a multiple of PBS_NCPUS - untested
CMDS_PATH=$HOME/tests/h2/getH2/PBS/cmds.txt

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file ${CMDS_PATH} --timeout 172800

# Combine output into a single file
cd /scratch/ht96/nb9894/h2/getH2

cat ./*_done.csv >> $SAVEDIR/out_h2.csv

# Delete loose files with
rm ./*_done.csv