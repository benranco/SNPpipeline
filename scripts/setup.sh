#!/bin/bash

echo "running setup.sh"

# (1) run the full pipeline
# (2) just process the data and do not generate reports (ie. just run the first half of the pipeline), 
# (3) just generate reports based on data that has already been processed by the first half of the pipeline (ie. just run the second half of the pipeline assuming the first half has already been run).
# (4) just generate reports beginning AFTER the filled_report.csv, assuming it has already been generated.
what_to_run=$1
generate_probability_report=$2

exitCode=0

if [[ what_to_run -ne 2 ]]
then
    Rscript ./scripts/install_R_packages.R $generate_probability_report
    exitCode=$?
fi

#make sure data isnt modified
chmod 555 ./data/* || echo "Just tried to change permissions of all the raw data to read+execute so as to avoid accidentally modifying it, but failed, which is fine because it means I don't have permissions to change their permissions, so I don't have permissions to actually modify their contents either."


exit $exitCode


