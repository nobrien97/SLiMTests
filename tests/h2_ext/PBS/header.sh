echo -e "gen,seed,modelindex,meanH,VA,phenomean,phenovar,dist,mean_w,deltaPheno,deltaW" | cat - out_stabsel_means.csv > d_means.csv

echo -e "gen,seed,modelindex,mutType,mutID,position,constraint,originGen,value,chi,Freq,mutCount,fixGen" | cat - out_stabsel_muts.csv > d_muts.csv

echo -e "gen,seed,modelindex,Q1_AUC,Q1_aZ,Q1_bZ,Q1_KZ,Q1_KXZ,Q2_AUC,Q2_aZ,Q2_bZ,Q2_KZ,Q2_KXZ,Q3_AUC,Q3_aZ,Q3_bZ,Q3_KZ,Q3_KXZ" | cat - out_stabsel_medodepar.csv > d_odepar.csv

echo -e "gen,seed,modelindex,vaZ,caZbZ,caZKZ,caZKXZ,vbZ,cbZKZ,cbZKXZ,vKZ,cKZKXZ,vKXZ" | cat - out_stabsel_GMat.csv > d_GMat.csv
echo -e "gen,seed,modelindex,Z_o,aZ_o,bZ_o,KZ_o,KXZ_o,Z_p,aZ_p,bZ_p,KZ_p,KXZ_p" | cat - out_stabsel_indPheno.csv > d_indPheno.csv

# Combine means and muts 
join -j1 -o1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12 <(<out_stabsel_muts.csv awk 'BEGIN { FS="," } ; {print $1"-"$2"-"$3" "$0}' | sort -k1,1) <(<out_stabsel_means.csv awk 'BEGIN { FS="," } ; {print $1"-"$2"-"$3" "$0}' | sort -k1,1) > d_combined.csv

awk 'NR==FNR{a[$1,$2,$3]=$1OFS$2;next}{$13=a[$1,$2,$3];print}' OFS=',' d_means.csv d_muts.csv > d_combined.csv
