#!/bin/bash -l

# This script removes the models we need to rerun from the done simulations folder in calcLD
# done files consist of ${MODELINDEX}_${SEED}
cd ~/tests/standingVar/calcLD

# First, get modelindices and seeds into an array

while IFS=, read -ra array; do
  MODELINDEX+=("${array[0]}")
  SEED+=("${array[1]}")
done < ./rerun_models.csv

# Then we need to iterate over the array and remove files with matching values
cd ./done
for (( i=0 ; i<${#MODELINDEX[@]} ; i++));
do
    rm ${MODELINDEX[i]}_${SEED[i]}
done