#!/bin/bash -l

module load python3/3.12.1

export PYTHONPATH=/g/data/ht96/nb9894/py_libs/lib/python3.12/site-packages:$PYTHONPATH

# Create a virtual environment if already not created.
python3 -m venv --system-site-packages /g/data/ht96/nb9894/py_libs/virtual/environment/ruggedness_3.12.1
 
# Activate your virtual environment
source /g/data/ht96/nb9894/py_libs/virtual/environment/ruggedness_3.12.1/bin/activate

python3 -m pip install -v --prefix /g/data/ht96/nb9894/py_libs --no-cache-dir pythae
python3 -m pip install -v --prefix /g/data/ht96/nb9894/py_libs --no-cache-dir pandas


deactivate