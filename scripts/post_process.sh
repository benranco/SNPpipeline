#!/bin/bash

input=$1 
single=$2
remove_indels=$3

export PERL5LIB=`pwd`"/tools/vcftools/src/perl"

echo "running post_process.sh"

filterInputFilename="$1.vcf"

if [[ $remove_indels == 1 ]]
then
  echo "removing indels"
  if [[ $single == 1 ]]
  then
      ./tools/vcftools/src/cpp/vcftools --vcf ./outputTemp/single/$1.vcf --remove-indels --recode --recode-INFO-all --out ./outputTemp/single/$1
  fi

  if [[ $single == 2 ]]
  then
      ./tools/vcftools/src/cpp/vcftools --vcf ./outputTemp/pooled/$1.vcf --remove-indels --recode --recode-INFO-all --out ./outputTemp/pooled/$1
  fi

  filterInputFilename="$1.recode.vcf"
fi

echo "filterInputFilename : "$filterInputFilename

echo "filtering using cutoffs"
if [[ $single == 1 ]]
then
	cat ./outputTemp/single/$filterInputFilename | ./tools/vcftools/src/perl/vcf-annotate --filter Qual=5 > ./outputTemp/single/$1_cutoff
fi

if [[ $single == 2 ]]
then
	cat ./outputTemp/pooled/$filterInputFilename | ./tools/vcftools/src/perl/vcf-annotate --filter Qual=5 > ./outputTemp/pooled/$1_cutoff
fi

echo "tab conversion"

if [[ $single == 1 ]]
then
    ./tools/vcftools/src/perl/vcf-to-tab < ./outputTemp/single/$1_cutoff > ./output/single/$1.tab
fi

if [[ $single == 2 ]]
then
    ./tools/vcftools/src/perl/vcf-to-tab < ./outputTemp/pooled/$1_cutoff > ./output/pooled/$1.tab
fi


