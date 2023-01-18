#!/bin/bash -l

# Setup files etc. for getH2_newTestCross.nsh

# This job uses a split version of slim_haplos for better speed
# To use, make sure you set up the right number of jobs for the number of input files (n-1, index starts at 0)
JOBNAME=getH2_newTestCross
JOBNAME_PREFIX=newTestCross/additive

# Make output folder
mkdir /scratch/user/uqnobri4/$JOBNAME_PREFIX/$JOBNAME
mkdir $HOME/tests/$JOBNAME_PREFIX/$JOBNAME/done

# Go to /g/data/ and split the input into 4 pieces - (11808*2000)/4 = 5904000: 
cd /scratch/user/uqnobri4/$JOBNAME_PREFIX
split -d -l 5904000 slim_haplo.csv 
