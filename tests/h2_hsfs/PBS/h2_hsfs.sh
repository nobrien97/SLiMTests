#!/bin/bash -l
#PBS -P ht96
#PBS -l walltime=24:00:00
#PBS -l ncpus=1440
#PBS -l mem=2000GB
#PBS -l jobfs=1000GB
#PBS -l storage=scratch/ht96+gdata/ht96
  
  
ECHO=/bin/echo

$ECHO "Creating outputs folders..."
# create our folders
cd $PBS_O_WORKDIR

# Make output folder
mkdir /scratch/ht96/nb9894/h2_hsfs
mkdir /g/data/ht96/nb9894/h2_hsfs
mkdir $HOME/tests/h2_hsfs/done

# Pre-job file manipulation goes here ...
# 
# INSERT CODE
cd $PBS_O_WORKDIR
SAVEDIR=/g/data/ht96/nb9894/h2_hsfs


# ========================================================================
# .... USER INSERTION OF EXECUTABLE LINE HERE 
# ========================================================================
# Make sure we're at the right place so we can find the bash script to run

$ECHO "Running nci-parallel..."
# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load nci-parallel/1.0.0a
export ncores_per_task=1
export ncores_per_numanode=12

# Calculate the range of parameter combinations we are exploring this job
# CAUTION: may error if CUR_TOT is not a multiple of PBS_NCPUS - untested
CMDS_PATH=$HOME/tests/h2_hsfs/PBS/cmds.txt

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file ${CMDS_PATH} --timeout 172800

$ECHO "All jobs finished, moving output..."

# Combine output into a single file
cd /scratch/ht96/nb9894/h2_hsfs/

cat ./slim_pos* >> $SAVEDIR/slim_pos.csv
cat ./slim_qg* >> $SAVEDIR/slim_qg.csv
cat ./slim_opt* >> $SAVEDIR/slim_opt.csv
cat ./slim_muts* >> $SAVEDIR/slim_muts.csv
cat ./slim_dict* >> $SAVEDIR/slim_dict.csv
cat ./slim_pedigree* >> $SAVEDIR/slim_pedigree.csv
cat ./slim_haplo* >> $SAVEDIR/slim_haplo.csv
cat ./slim_sampled_pheno* >> $SAVEDIR/slim_sampled_pheno.csv
cat ./slim_genmap* >> $SAVEDIR/slim_genmap.csv
cat ./slim_fx* >> $SAVEDIR/slim_fx.csv
cat ./slim_time* >> $SAVEDIR/slim_time.csv

# Zip LD matrices
/bin/zip -q -Z bzip2 $SAVEDIR/slim_ld.zip ./slim_ld* 

# Delete loose files with seed and model indices
find -regex ".*[0-9]*_*[0-9].csv+" -delete
rm *.tsv
