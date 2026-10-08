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


# Start h2o java instance 
ip=$(hostname -I | awk '{print $1}')
java -Xmx40g -jar $HOME/R/x86_64-pc-linux-gnu-library/4.0/h2o/java/h2o.jar -ip $ip -port 12345 -quiet > /dev/null 2>&1 &
h2oid=$!
sleep 15 # Sleep to let h2o start up

# Calculate landscape
module load R/4.0.0
RSCRIPTNAME=$TESTDIR/R/pca_landscape.R
Rscript ${RSCRIPTNAME} ${ip}

# Close h2o
kill $h2oid
