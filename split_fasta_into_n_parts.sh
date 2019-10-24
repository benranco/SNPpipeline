#!/bin/bash

#####################################################
# This is a helper script to split a fasta file into a given number of equal 
# pieces (according to number of sequences in each piece).
# 
# To run this script, edit the below input parameters, then from a command-line
# Terminal execute either:
#      bash ./split_fasta_into_n_parts.sh
# or simply:
#     ./split_fasta_into_n_parts.sh
#
#####################################################
# Input parameters:

# The path to the folder containing fasta reference file you want to split
path="/path/to/folder"

# The name of the fasta reference file you want to split
inputFastaFile="filename.fasta"

# The number of parts you want to split the fasta reference file into
numParts=10

# The first part of the file name for the output files
outFileNamePrefix="filename"

#####################################################
# Execution code (nothing below here needs to be edited):

totalSequences=$(grep -c "^>" $path/$inputFastaFile )
nSequences=$(expr $totalSequences / $numParts + 1)

echo "-------------------"
echo "The total number of sequences in "$inputFastaFile" is "$totalSequences"."
echo "Splitting the file into "$numParts" files of "$nSequences" sequences each (mabe fewer in the last file)..."

# This command does all the work. The -v indicates an input variable. The stuff inside 
# the '' is the code to execute.
awk -v numSequences=$nSequences -v outFName=$path/$outFileNamePrefix -v postfix=fasta '
   /^>/ { n++; if (n % numSequences == 1) { i++; close(fname); fname = sprintf("%s.%02d.%08d.%s" , outFName, i, n, postfix) } }
   { print >> fname }
' $path/$inputFastaFile

echo "done."
