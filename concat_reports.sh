#!/bin/bash

#####################################################
# This is a helper script designed to concatenate .csv reports sharing an identical name and 
# identical header lines in multiple folders into one report. This script will 
# operate on all .csv files in the given input folders, and will assume that files with matching 
# names are intended to be concatenated together.
#
# This script is intended for those situations when the reference fasta file used for the SNPpipeline
# is split up into pieces, and the pipeline is run concurrently on each of these pieces. The final reports
# generated from each of these concurrent pipelines can then be concatenated together.
#
# To run this script, modify the input parameters below before running the script.

#####################################################
# Input Parameters: modify these parameters before running the script

# The root path. All subfolders are assumed to be relative to this path.
path="."

# The subfolder containing the reports that are to be used as the masters, which will be concatenated onto.
# If you wish to keep an unmodified version of the originals, make sure you save copies first.
# We assume that all the .csv files in the master folder are reports, and that all of them are matched by 
# equivalently named reports in other folders which will be concatenated onto these master reports.
# only .csv files will be operated on. All other files will be ignored.
masterFolder="master"

# An array of all subfolders from which to draw reports that will be concatenated onto the master reports.
inputFolders=(reports02 reports03 reports04 reports05 reports06 reports07 reports08 reports09 reports10)


####################################################
# Execution:

echo "Running concat_reports.sh"
echo "Reports to concatenate: "
ls $path/$masterFolder | grep ".csv"

for reportName in $(ls $path/$masterFolder | grep ".csv")
do
  masterReport=$path/$masterFolder/$reportName

  echo " "
  echo "============================================="
  echo "Concatenating "$masterReport
  masterLines=`wc -l $masterReport`
  echo "masterReport num lines: "$masterLines
  echo "All reports concatenated to this report will be concatenated without their header line."
  echo "------------------------------"

  for folder in ${inputFolders[@]}
  do
    subReport=$path/$folder/$reportName
    
    cmp <(head -n 1 $masterReport) <(head -n 1 $subReport)
    # the $? stores the exit code of the last command, which should be zero if
    # the first line of each of the files in the above cmp command (which is like diff) are identical.
    if [[ $? == 0 ]] 
    then
      lines=`wc -l $subReport`
      echo "  Adding: (num lines, incl header = "$lines")"
      tail -n +2 $subReport >> $masterReport
      masterLines=`wc -l $masterReport`
      echo "  masterReport num lines: "$masterLines
    else
      echo "!!! Did not concatenate "$subReport" to "$masterReport" because the headers didn't match or one of the files doesn't exist."
    fi
  done # end inner for-loop

done # end outer for-loop

echo " "
echo "============================================="
echo "Finished concat_reports.sh"

