# This is intended to pre-install all packages required by the R scripts in the SNPpipeline. 
# The idea is to call this script at the beginning of the pipeline so if there are any 
# problems in installing them we find out right away.
#
# In order to successfully install these you might need to run R with root permissions, 
# depending on how R was installed.

args <- commandArgs(trailingOnly = TRUE)
generate_probability_report <- args[1]
generate_probability_report <- as.integer(generate_probability_report)

message("****************************************************")
message("Installing R packages used in the pipeline.")
message(" ")

# Used in report_gen_p2.R
if(!require(seqinr)) {
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}
library(seqinr)

# Used in getDepthStats-parallel.R
if(!require(doParallel)) {
  install.packages('doParallel', repos='http://cran.us.r-project.org')
}  
library(doParallel)  

if(generate_probability_report) {
  # Used in report_gen_p2.R
  # This install proceedure is tested to work in R 4.2.1, specifically: BiocManager::install(version = "3.15"
  # It may require keyboard input in response to a question (select "a" for all).
  # If there are errors (you can ignore warnings) during installation, it might be because some required 
  # linux libraries need to be installed first. You can figure that out be reading through the error messages.
  # For example, I needed to install libxml2-devel and openssl-devel via dnf.  
  if(!require(VariantAnnotation)) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
    }
    BiocManager::install(version = "3.15")
    BiocManager::install("VariantAnnotation")
  }
  library(VariantAnnotation)
}


message("Finished installing R packages used in the pipeline.")
message("****************************************************")
message(" ")

