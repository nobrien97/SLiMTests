#!/bin/bash -l

# Splits mutations file to limit epistasis calculation memory usage
JOB_PATH=/mnt/c/GitHub/SLiMTests/tests/standingVar/analysis/R/
DATA_PATH=/mnt/d/SLiMTests/tests/standingVar/dec11/
FILENAME=slim_muts_sbst.csv
DATA=${DATA_PATH}/${FILENAME}

# Remove lines with gens < 50000
awk -F, '$1 >= 49500' $DATA > slim_muts_sbst_sml.csv


# Split the data according to the seed and modelindex
awk -F, '{print > ("muts_" $2"_" $3 ".csv")}' $DATA

function print_buffer {
  # logic to extract $2 and $3 from the buffer and create the filename
  print "$1" > "muts_" $2"_" $3 ".csv"
}
awk -F, '{ buffer = buffer $0 ORS } NR % 1000 == 0 { print_buffer buffer $2 $3; buffer="" } END { print_buffer(buffer) }' $DATA
