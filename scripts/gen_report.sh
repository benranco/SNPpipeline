#!/bin/bash
echo "running gen_report.sh"

mafcutoff=$1
generate_chi_sq_report=$2 
generate_probability_report=$3 
generate_depth_stats_report=$4
haploidOrDiploid=$5

Rscript ./scripts/report_gen.R `pwd`
Rscript ./scripts/report_subset.R `pwd`

counter=0
ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -4)*2 ))"

echo "running parallel processing"

for f in $(ls ./reporttemp | grep -v "filled.Rds")
do    
	if [ $counter -ge $ncore ]; then	
		wait
		#echo "wait reset"
		counter=0
	fi

	Rscript ./scripts/parallel_process.R `pwd` $f $haploidOrDiploid &
	counter=$(( $counter + 1 ))
done
wait

Rscript ./scripts/report_gen_p2.R `pwd` $mafcutoff $generate_chi_sq_report $generate_probability_report $haploidOrDiploid

if [[ generate_depth_stats_report -eq 1 ]]
then
    Rscript ./scripts/getDepthStats-parallel.R `pwd`
fi


