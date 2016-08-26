# SNPpipeline
discovery of SNPs

pipeline developed with bowtie2 for alignment, freebayes for SNP detection, hard coded variables and various connecting tools 

Developed for Dr.Jun-Jun Liu and Dr.Isabel Leal at the Pacific Forestry Center, NRCAN, Government of Canada

## Requirements

Install dependencies using the following commands (requires admin privledges)

`sudo yum update`

`sudo yum install epel-release`

`sudo yum install R libcurl-devel.x86_64 libxml2-devel.x86_64` 

Install R dependecies should be automatic through the script, although some screens might require user input. The answer is always update all, and use personal library

## How to run pipeline 

only run this pipeline with copies of data, dont use one copy of data for everything

put reference, or references into the reference folder.
if there are multiple reference files they will be merged into one file
also dont name your reference formatted_reference.fasta, because if it is it'll assume that the reference step has already been run

put samples into data. they must be labeled as such or the pipeline wont work

<some identifier>_R1.<some fasta format>
<some identifier>_R2.<some fasta format>

examples of valid fasta formats are .fa, .fasta, .mfa, .fna
compressed fasta files such as fastq.gz also work

make sure all reads are paired, otherwise the alignment will not work 

edit the start.sh with correct arguments in start.sh, and then run start.sh 

report generation:

`report.csv` is the initial report from all SNPs merged from all samples

`filled_report.csv` is the report with NAs filled in with either reference, or NA for no alignment 

`edited_report.csv` is report removing SNP sites with mostly NA, or 1 genotype 

`percentage_snps.csv` is where MAF is calculated, without MAF cutoff 

`mutation_percentage.csv` is percentage of snps divided by length of site, using sites making the MAF cutoff 

`MAF_cutoff_report.csv` is `edited_report.csv` with sites meeting the MAF cutoff

`MAF_cutoff_report_chi.csv` is `MAF_cutoff_report.csv` with the genotype with most observations replaced with "H" and the second most genotype with "A"

`probabilities.csv` is `MAF_cutoff_report.csv`  with extracted percentages from vcf files
