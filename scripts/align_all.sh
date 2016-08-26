#!/bin/bash
#align all files

echo "running align_all.sh"

input=$1

ref="./referenceTemp/formatted_output"


if [[ $3 == 0 ]]; then
    
    echo "making comma seperated list"
    data=""
    for one in $(ls ./data)
    do
        if [[ -z $data ]]; then
	    data="./data/$one"
        else
	    data+=",./data/$one"
        fi
    done

    echo $data

    #align with bowtie2, save to temp
    echo "alignment with bowtie" 
	./tools/bowtie2-2.2.9/bowtie2 -p $2 -x $ref -U $data --local --very-sensitive-local -S "./dataTemp/pooled/$input.sam"

elif [[ $3 == 1 ]]; then

    echo "making comma seperated lists for R1 and R2"
    dataone=""
    for one in $(ls ./data | grep _R1)
    do
        if [[ -z $dataone ]]; then
	    dataone="./data/$one"
        else
	    dataone+=",./data/$one"
        fi
    done

    datatwo=""
    for two in $(ls ./data | grep _R2)
    do
        if [[ -z $datatwo ]]; then
	    datatwo="./data/$two"
        else
	    datatwo+=",./data/$two"
        fi
    done

    echo $dataone
    echo $datatwo

    #align with bowtie2, save to temp
    echo "alignment with bowtie" 
    ./tools/bowtie2-2.2.9/bowtie2 -p $2 -x $ref -1 $dataone -2 $datatwo --local --very-sensitive-local -S "./dataTemp/pooled/$input.sam"
fi

echo "--------------------"
