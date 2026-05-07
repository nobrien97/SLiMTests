#!/bin/bash -l

# Only run r = 0.1 models, write done files for the others
# the model indices we want to run are 11-15
# So we need done files for 1-10

NJOB_COUNTER=0
NJOB=0
NJOB_FILES=1560
for i in {1..10}
do
    for j in {1..208}
    do
        touch ${i}_${j}_${NJOB}
        ((NJOB_COUNTER+=1))

        if (( $NJOB_COUNTER >= $NJOB_FILES )); then
            NJOB=1
        fi
    done
done



