#!/bin/bash -l

cd /g/data/ht96/nb9894/newMotifs/paper1/randomisedStartsM

# Combine data
cat slim_pos_0.csv slim_pos_1.csv >> slim_pos.csv
cat slim_opt_0.csv slim_opt_1.csv >> slim_opt.csv
cat slim_muts_0.csv slim_muts_1.csv >> slim_muts.csv
cat slim_qg_0.csv slim_qg_1.csv >> slim_qg.csv
cat slim_haplo_fix_0.csv slim_haplo_fix_1.csv >> slim_haplo_fix.csv
cat slim_haplo_0.csv slim_haplo_1.csv >> slim_haplo.csv
cat slim_sampled_pheno_0.csv slim_sampled_pheno_1.csv >> slim_sampled_pheno.csv
cat slim_indPheno_0.csv slim_indPheno_1.csv >> slim_indPheno.csv
cat slim_sampled_moltrait_0.csv slim_sampled_moltrait_1.csv >> slim_sampled_moltrait.csv
cat slim_fx_0.csv slim_fx_1.csv >> slim_fx.csv
cat slim_locusHo_0.csv slim_locusHo_1.csv >> slim_locusHo.csv
cat slim_PMmat_0.csv slim_PMmat_1.csv >> slim_PMmat.csv
cat slim_relPos_0.csv slim_relPos_1.csv >> slim_relPos.csv
cat slim_relVals_0.csv slim_relVals_1.csv >> slim_relVals.csv
cat slim_sharedmutfreqs_0.csv slim_sharedmutfreqs_1.csv >> slim_sharedmutfreqs.csv
cat slim_mutvar_0.csv slim_mutvar_1.csv >> slim_mutvar.csv
cat slim_mutvar_percomp_0.csv slim_mutvar_percomp_1.csv >> slim_mutvar_percomp.csv

# Remove split datasets
# find -regex ".*[0-1].csv" -delete
