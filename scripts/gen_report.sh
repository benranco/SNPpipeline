#!/bin/bash
echo "running gen_report.sh ---" `date`

mafcutoff=$1
generate_chi_sq_report=$2 
generate_probability_report=$3 
generate_depth_stats_report=$4
haploidOrDiploid=$5
wereIndelsRemoved=$6
# (1) run the full pipeline
# (2) just process the data and do not generate reports (ie. just run the first half of the pipeline), 
# (3) just generate reports based on data that has already been processed by the first half of the pipeline (ie. just run the second half of the pipeline assuming the first half has already been run).
# (4) just generate reports beginning AFTER the filled_report.csv, assuming it has already been generated.
what_to_run=$7

if [[ what_to_run -ne 4 ]]
then
  echo "generating report.csv"
  Rscript ./scripts/report_gen.R `pwd`
  Rscript ./scripts/report_subset.R `pwd`

  counter=0
  ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -4)*2 ))"

  echo "running parallel processing to generate filled_report.csv"

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
fi

Rscript ./scripts/report_gen_p2.R `pwd` $mafcutoff $generate_chi_sq_report $generate_probability_report $haploidOrDiploid $wereIndelsRemoved $what_to_run


if [[ generate_depth_stats_report -eq 1 ]]
then
    echo "getting depth statistics ---" `date`
    Rscript ./scripts/getDepthStats-parallel.R `pwd`
fi

echo "finished gen_report.sh ---" `date`

