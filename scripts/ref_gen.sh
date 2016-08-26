#!/bin/bash

ncore=$1
export PERL5LIB=$PERL5LIB:`pwd`/tools/vcftools/src/perl

echo "running ref_gen.sh"

rm -f ./reference/*.fai

c=0
for f in $(ls "./reference")
do 
   ((c++))
   ref="$f"
done 

if [[ $c > 1 ]]
then
    echo "found multiple references, concat to one file"
    cat ./reference/* >> ./reference/output.fasta
    for df in $(ls ./reference | grep -v output.fasta) 
    do
	rm ./reference/$df
    done
    ref="output.fasta"    
elif [[ $c == 0 ]]
then
    echo "found no reference"
    exit 1
else
    echo "found one reference"
fi

if [ "$ref" != "formatted_output.fasta" ]
then
    echo "formatting fasta file to consistent line lengths"
    java -jar ./tools/picard-tools-2.3.0/picard.jar NormalizeFasta I=./reference/$ref O=./reference/"formatted_output.fasta" LINE_LENGTH=60 
    #mv ./reference/$ref ./reference/"formatted_output.fasta"
fi

echo "remove original fasta file"
for dfa in $(ls ./reference | grep -v "formatted_output.fasta")
do
    rm ./reference/$dfa
done

ref="./referenceTemp/formatted_output"
newref="./reference/formatted_output.fasta"

echo "build new ref"
./tools/bowtie2-2.2.9/bowtie2-build --threads $ncore $newref $ref 

echo "build fai for reference"
./tools/samtools-1.3.1/samtools faidx "./reference/formatted_output.fasta"
