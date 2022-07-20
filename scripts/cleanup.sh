#!/bin/bash

echo "running cleanup.sh"

# (1) run the full pipeline
# (2) just process the data and do not generate reports (ie. just run the first half of the pipeline), 
# (3) just generate reports based on data that has already been processed by the first half of the pipeline (ie. just run the second half of the pipeline assuming the first half has already been run).
# (4) just generate reports beginning AFTER the filled_report.csv, assuming it has already been generated.
what_to_run=$7

if [[ what_to_run -eq 1 ]] || [[ what_to_run -eq 2 ]]
then
  rm -f ./dataTemp/single/*  
  rm -f ./dataTemp/pooled/*  
  rm -f ./dataTemp/picardTmp/*
  rm -f ./dataTemp/depthfiles/*
  rm -f ./output/single/*  
  rm -f ./output/pooled/*   
  rm -f ./outputTemp/single/*  
  rm -f ./outputTemp/pooled/*  
  rm -f ./logs/single/*  
  rm -f ./logs/pooled/*   
  rm -f ./logs/*
  rm -f ./output/*
  rm -f ./referenceTemp/*
fi

if [[ what_to_run -eq 1 ]] || [[ what_to_run -eq 3 ]]
then
  rm -f ./reporttemp/*  
  rm -f ./reports/*.csv  
  rm -f ./reports/*.linkage  
fi
 
if [[ what_to_run -eq 4 ]]
then

  mv ./reports/report.csv ./reports/report.keep
  mv ./reports/filled_report.csv ./reports/filled_report.keep
    
  rm -f ./reports/*.csv  
  rm -f ./reports/*.linkage  
  
  mv ./reports/report.keep ./reports/report.csv
  mv ./reports/filled_report.keep ./reports/filled_report.csv

fi


