#!/bin/bash
echo "running gen_report.sh"

Rscript ./scripts/report_gen.R `pwd`
Rscript ./scripts/report_subset.R `pwd`

counter=0
ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -4)*2 ))"

echo "running parallel processing"

for f in $(ls ./reporttemp | grep -v "filled.Rds")
do    
	if [ $counter -lt $ncore ]; then	
		Rscript ./scripts/parallel_process.R `pwd` $f &
		counter=$(( $counter + 1 ))
	else
		wait
		#echo "wait reset"
		counter=0
	fi
done
wait

Rscript ./scripts/report_gen_p2.R `pwd` $1
