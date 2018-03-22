#!/bin/bash

echo "running setup.sh"

# what_to_run: Parameter specified in start.sh: (1) run the full pipeline, (2) just process the data and do not generate reports (ie. just run the first half of the pipeline), (3) just generate reports based on data that has already been processed by the first half of the pipeline (ie. just run the second half of the pipeline assuming the first half has already been run).
what_to_run=$1

exitCode=0

if [[ what_to_run -eq 1 ]] || [[ what_to_run -eq 3 ]]
then
    sudo Rscript ./scripts/install_R_packages.R
    exitCode=$?
fi

exit $exitCode

