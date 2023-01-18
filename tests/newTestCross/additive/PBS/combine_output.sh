#!/bin/bash
JOBNAME=newTestCross


# Combine output into a single file
cd /scratch/user/uqnobri4/$JOBNAME/additive/

cat ./slim_pos* >> ./slim_pos.csv
cat ./slim_qg* >> ./slim_qg.csv
cat ./slim_opt* >> ./slim_opt.csv
cat ./slim_muts* >> ./slim_muts.csv
cat ./slim_dict* >> ./slim_dict.csv
cat ./slim_haplo* >> ./slim_haplo.csv
cat ./slim_sampled_pheno* >> ./slim_sampled_pheno.csv
cat ./slim_genmap* >> ./slim_genmap.csv
cat ./slim_fx* >> ./slim_fx.csv
cat ./slim_time* >> ./slim_time.csv

# Zip LD matrices
/bin/zip -q -Z bzip2 ./slim_ld.zip ./slim_ld* 

# Delete loose files with seed and model indices
find -regex ".*[0-9]*_*[0-9].csv+" -delete
rm *.tsv
