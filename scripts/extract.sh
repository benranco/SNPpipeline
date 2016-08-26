#!/bin/bash
#get output from temp, find SNPs using FreeBayes

echo "running extract.sh"

if [ $3 == 1 ]
then
    path="single"
elif [ $3 == 2 ]
then
    path="pooled"
fi

input="./dataTemp/$path/$1.sam"
 
echo "generate bam"
./tools/samtools-1.3.1/samtools view -Sb $input > ./dataTemp/$path/$1.bam   
rm -f $input 

echo "sort bam"
./tools/samtools-1.3.1/samtools sort -o ./dataTemp/$path/$1_sorted.bam ./dataTemp/$path/$1.bam 
rm -f ./dataTemp/$path/$1.bam

echo "index bam"
./tools/samtools-1.3.1/samtools index ./dataTemp/$path/$1_sorted.bam 

binput="./dataTemp/$path/$1_sorted.bam"
#echo $binput
ref="./reference/formatted_output.fasta"
#echo $1

echo "find SNPs with freebayes"

if [ $4 == 0 ]
then
	./tools/freebayes/bin/freebayes -f $ref -p $2 --pvar 0.75 --min-alternate-count 1 --theta 0.01  --use-best-n-alleles 4 --min-alternate-fraction 0.8 $binput 1> ./outputTemp/$path/$1.vcf
elif [ $4 == 1 ]
then 
	./tools/freebayes/bin/freebayes -f $ref -p $2 $binput 1> ./outputTemp/$path/$1.vcf
fi
