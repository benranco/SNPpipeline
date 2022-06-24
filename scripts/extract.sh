#!/bin/bash
#get output from temp, find SNPs using FreeBayes

echo "running extract.sh ---" `date`

# the $1 input parameter is the input filename minus the extension
ploidy=$2
single=$3
default=$4

if [ $single == 1 ]
then
    path="single"
elif [ $single == 2 ]
then
    path="pooled"
fi

input="./dataTemp/$path/$1.sam"
 
echo "generate bam ---" `date`
./tools/samtools-1.3.1/samtools view -Sb $input > ./dataTemp/$path/$1.bam   
rm -f $input 

echo "sort bam ---" `date`
./tools/samtools-1.3.1/samtools sort -o ./dataTemp/$path/$1_sorted.bam ./dataTemp/$path/$1.bam 
rm -f ./dataTemp/$path/$1.bam


echo "mark duplicates ---" `date`

maxf=$(( ($(ulimit -n) - 24) ))
if [ -z "$maxf" ]; then maxf=1000; fi

java -Djava.io.tmpdir=./dataTemp/picardTmp -jar ./tools/picard-2.20.7/picard.jar MarkDuplicates TMP_DIR=./dataTemp/picardTmp  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$maxf I=./dataTemp/$path/$1_sorted.bam O=./dataTemp/$path/$1_sorted_markDup.bam M=./dataTemp/$path/$1_sorted_markDup_metrics.txt > ./dataTemp/$path/$1_sorted_markDup_log.txt
rm -f ./dataTemp/$path/$1_sorted.bam

echo "index bam ---" `date`
./tools/samtools-1.3.1/samtools index ./dataTemp/$path/$1_sorted_markDup.bam 

binput="./dataTemp/$path/$1_sorted_markDup.bam"

#echo $binput
ref="./reference/formatted_output.fasta"
#echo $1

echo "find SNPs with freebayes ---" `date`

if [ $default == 0 ]
then
  # parameters from a research paper
	./tools/freebayes/bin/freebayes -f $ref -p $ploidy --pvar 0.75 --min-alternate-count 1 --theta 0.01  --use-best-n-alleles 4 --min-alternate-fraction 0.8 $binput 1> ./outputTemp/$path/$1.vcf
elif [ $default == 1 ]
then 
  # default params for single raw data files
	./tools/freebayes/bin/freebayes -f $ref -p $ploidy $binput 1> ./outputTemp/$path/$1.vcf
elif [ $default == 2 ]
then 
  # default params for pooled raw data files
	./tools/freebayes/bin/freebayes -f $ref -p $ploidy --min-alternate-fraction 0.001 $binput 1> ./outputTemp/$path/$1.vcf
fi

echo "finished ---" `date`
