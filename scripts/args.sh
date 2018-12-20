#!/bin/bash
#get arguments here to call other scripts

datapoint=$1
name=$2
single=$3
ploidy=$4
ncore=$5
default=$6
paired=$7
remove_indels=$8

#calling methods required for future calls 

#make sure data isnt modified
chmod 555 ./data/*

if [ $single == 1 ]
then
    bash ./scripts/align.sh $datapoint $ncore $paired > ./logs/single/$datapoint"_align.log" 2>&1 
elif [ $single == 2 ]
then
    bash ./scripts/align_all.sh $datapoint $ncore $paired > ./logs/pooled/$datapoint"_align.log" 2>&1
fi

if [ $single == 1 ]
then
    bash ./scripts/extract.sh $datapoint $ploidy $single $default > ./logs/single/$datapoint"_extract.log" 2>&1
elif [ $single == 2 ]
then
    bash ./scripts/extract.sh $datapoint $ploidy $single $default > ./logs/pooled/$datapoint"_extract.log" 2>&1
fi

if [ $single == 1 ]
then
    bash ./scripts/post_process.sh $datapoint $single $remove_indels > ./logs/single/$datapoint"_post_process.log" 2>&1
elif [ $single == 2 ]
then
    bash ./scripts/post_process.sh $datapoint $single $remove_indels > ./logs/pooled/$datapoint"_post_process.log" 2>&1
fi

