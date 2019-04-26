#!/bin/bash
echo "running gen_report.sh ---" `date`

mafcutoff=$1
generate_chi_sq_report=$2 
generate_probability_report=$3 
generate_depth_stats_report=$4
haploidOrDiploid=$5
wereIndelsRemoved=$6

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

echo "finished parallel processing ---" `date`

Rscript ./scripts/report_gen_p2.R `pwd` $mafcutoff $generate_chi_sq_report $generate_probability_report $haploidOrDiploid $wereIndelsRemoved


if [[ generate_depth_stats_report -eq 1 ]]
then
    echo "getting depth statistics ---" `date`
    Rscript ./scripts/getDepthStats-parallel.R `pwd`
fi

echo "finished gen_report.sh ---" `date`

