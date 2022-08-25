#!/bin/bash

# For further documentation on this pipeline, see the README.md file.

#call this pipeline something unique, replace spaces with underlines
name="junjun_pipeline"

# paired input raw data files (e.g. R1 & R2) vs unpaired files: (1) paired, (0) unpaired
paired=1

# whether each sample has its own raw data file(s), or samples are pooled into combined input file(s): (1) single, (2) pooled, (3) both
single=1

#ploidy values for single and pooled runs: (1) haploid, (2) diploid
singleploidy=1
pooledploidy=2


#use default freebayes parameters or not: (3) custom (not implemented), (1) default params for single raw data files, (2) default params for pooled raw data files (only difference is --min-alternate-fraction is set to 0.001 instead of the default of 0.2, (0) parameters from research paper
freebayesParams=1

#maf cutoff 
mafcut=0.05

# (1) remove indels after finding SNPs, (0) don't remove indels after finding SNPs
remove_indels=0

# (1) run the full pipeline
# (2) just process the data and do not generate reports (ie. just run the first half of the pipeline), 
# (3) just generate reports based on data that has already been processed by the first half of the pipeline (ie. just run the second half of the pipeline assuming the first half has already been run).
# (4) just generate reports beginning AFTER the filled_report.csv, assuming it has already been generated.
what_to_run=1

# some optional reports which you might choose not to generate in order to save time:

# (1) yes, (0) no
generate_chi_sq_report=1

# (1) yes, (0) no. This can be very time consuming, and as yet has not been useful.
generate_probability_report=0

# (1) yes, (0) no. This can be very time consuming.
generate_depth_stats_report=1

# The pipeline will assume it can use all available cpu's for the bowtie alignment except the number you specify here (to make it easier to coexist with other processes you might be running on the computer).
freeCpus=1

# Scroll down to see the code.
























########################################################################
# function formatPairedFileNames()
# This gets the list of the R1 and R2 raw data input file names, chops  
# off the end portion according to the various allowed input file types, 
# and returns in the $datapoints variable the unique id portion of the
# file names, after eliminating adjacent duplicates due to the R1 and
# R2 files sharing the same base name.
# Functions in bash are a little crazy, so the input parameter needs to 
# be $datapoints. So you'd call it like:
# datapoints=""
# formatFileNames $datapoints
formatPairedFileNames()
{
datapoints=$(ls "./data" \
              | sed "s/_R1.fastq.gz//g" \
              | sed "s/_R2.fastq.gz//g" \
              | sed "s/_R1.fasta.gz//g" \
              | sed "s/_R2.fasta.gz//g" \
              | sed "s/_R1.mfa.gz//g" \
              | sed "s/_R2.mfa.gz//g" \
              | sed "s/_R1.fna.gz//g" \
              | sed "s/_R2.fna.gz//g" \
              | sed "s/_R1.fa.gz//g" \
              | sed "s/_R2.fa.gz//g" \
              | sed "s/_R1.fq.gz//g" \
              | sed "s/_R2.fq.gz//g" \
              | sed "s/_R1.fastq//g" \
              | sed "s/_R2.fastq//g" \
              | sed "s/_R1.fasta//g" \
              | sed "s/_R2.fasta//g" \
              | sed "s/_R1.mfa//g" \
              | sed "s/_R2.mfa//g" \
              | sed "s/_R1.fna//g" \
              | sed "s/_R2.fna//g" \
              | sed "s/_R1.fa//g" \
              | sed "s/_R2.fa//g" \
              | sed "s/_R1.fq//g" \
              | sed "s/_R2.fq//g" \
              | uniq)
}


./scripts/setup.sh $what_to_run $generate_probability_report
# if setup had a bad exit code, abort the pipeline
if [ $? -gt 0 ]
then
    echo "Failed in setup. Aborting pipeline."
    exit $?
fi

./scripts/cleanup.sh $what_to_run

if [ $single -eq 3 ]; then
	ncore="$(( ($(grep -c ^processor /proc/cpuinfo) - $freeCpus)/2 ))"
else
	ncore="$(( ($(grep -c ^processor /proc/cpuinfo) - $freeCpus)))"
fi

if [ $ncore -lt 1 ]
then
    ncore=1
fi

# running the first half of the pipeline
if [[ what_to_run -eq 1 ]] || [[ what_to_run -eq 2 ]]
then
    ./scripts/ref_gen.sh $ncore > ./logs/$name"_reference_generation.log" 2>&1

    if [[ single -eq 1 ]]
    then
        echo "running single"
        if [[ paired -eq 0 ]]; then
            for datapoint in $(ls "./data"); do
                echo "Processing" $datapoint "---" `date`
                ./scripts/args.sh $datapoint $name $single $singleploidy $ncore $freebayesParams $paired $remove_indels
            done
        elif [[ paired -eq 1 ]]; then
            datapoints=""
            formatPairedFileNames $datapoints
            for datapoint in $datapoints
            do
                echo "Processing" $datapoint "---" `date`
                ./scripts/args.sh $datapoint $name $single $singleploidy $ncore $freebayesParams $paired $remove_indels
                #sync
                #echo 1 > /proc/sys/vm/drop_caches
            done
        fi
    elif [[ single -eq 2 ]]
    then
        echo "running pooled"
        ./scripts/args.sh $name $name $single $pooledploidy $ncore $freebayesParams $paired $remove_indels
    elif [[ single -eq 3 ]]
    then
        echo "running pooled in background"
        ./scripts/args.sh $name $name 2 $pooledploidy $ncore $freebayesParams $paired $remove_indels &
        echo "running single"
        if [[ paired -eq 0 ]]; then
            for datapoint in $(ls "./data"); do
                echo "Processing" $datapoint "---" `date`
                ./scripts/args.sh $datapoint $name 1 $singleploidy $ncore $freebayesParams $paired $remove_indels 
            done
        elif [[ paired -eq 1 ]]; then
            datapoints=""
            formatPairedFileNames $datapoints
            for datapoint in $datapoints
            do
                echo "Processing" $datapoint "---" `date`
                ./scripts/args.sh $datapoint $name 1 $singleploidy $ncore $freebayesParams $paired $remove_indels
                #sync
                #echo 1 > /proc/sys/vm/drop_caches
            done
        fi
        wait
    fi

    echo "Finished processing. Have not yet generated reports (scripts/gen_report.sh). ---" `date`
fi

# running the second half of the pipeline, report generation
if [[ what_to_run -eq 1 ]] || [[ what_to_run -eq 3 ]] || [[ what_to_run -eq 4 ]]
then
    echo "Generating reports. ---" `date`
    ./scripts/gen_report.sh $mafcut $generate_chi_sq_report $generate_probability_report $generate_depth_stats_report $singleploidy $remove_indels $what_to_run
fi


echo  -e "\033[33;5;7mPipeline $name has finished\033[0m"

date

 

