#!/bin/bash
#call this pipeline something unique, replace spaces with underlines
name="junjun_fluidigm"

#(1) single, (2) pooled, (3) both
single=1

#contact info for pipeline completion, doesnt work
email=""

#default freebayes parameters or not: (3) custom (not implemented), (1) default, (0) paper
default=1

#maf cutoff 
mafcut=0.3

# paired files (e.g. R1 & R2) vs unpaired files: (1) paired, (0) unpaired
paired=1

#ploidy values for single and pooled runs: (1) haploid, (2) diploid
singlep=1
pooledp=2


# Scroll down to see the code.































if [ $single -eq 3 ]; then
	ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -4)/2 ))"
else
	ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -4)))"
fi

if [ $ncore -lt 1 ]
then
    ncore=1
fi

./scripts/cleanup.sh
./scripts/ref_gen.sh $ncore > ./logs/$name"_reference_generation.log" 2>&1

if [[ single -eq 1 ]]
then
    echo "running single"
    if [[ paired -eq 0 ]]; then
        for datapoint in $(ls "./data"); do
            ./scripts/args.sh $datapoint $name $single $singlep $ncore $default $paired
        done
    elif [[ paired -eq 1 ]]; then
        for datapoint in $(ls "./data" | rev | cut -c 13- | rev | uniq)
        do
            ./scripts/args.sh $datapoint $name $single $singlep $ncore $default $paired
            #sync
            #echo 1 > /proc/sys/vm/drop_caches
        done
    fi
elif [[ single -eq 2 ]]
then
    echo "running pooled"
    ./scripts/args.sh $name $name $single $pooledp $ncore $default $paired
elif [[ single -eq 3 ]]
then
    echo "running pooled in background"
    ./scripts/args.sh $name $name 2 $pooledp $ncore $default $paired &
    echo "running single"
    if [[ paired -eq 0 ]]; then
        for datapoint in $(ls "./data"); do
            ./scripts/args.sh $datapoint $name 1 $singlep $ncore $default $paired 
        done
    elif [[ paired -eq 1 ]]; then
        for datapoint in $(ls "./data" | rev | cut -c 13- | rev | uniq)
        do
            ./scripts/args.sh $datapoint $name 1 $singlep $ncore $default $paired
            #sync
            #echo 1 > /proc/sys/vm/drop_caches
        done
    fi
    wait
fi

echo "Finished processing. Have not yet generated reports (scripts/gen_report.sh)."

echo "Generating reports."
./scripts/gen_report.sh $mafcut

echo  -e "\033[33;5;7mPipeline $name has finished\033[0m"

date

 

