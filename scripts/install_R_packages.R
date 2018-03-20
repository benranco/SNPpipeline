# This is intended to pre-install all packages required by the R scripts in the SNPpipeline. 
# The idea is to call this script at the beginning of the pipeline so if there are any 
# problems in installing them we find out right away.
#
# In order to successfully install these you might need to run R with root permissions.

message("****************************************************")
message("Installing R packages used in the pipeline.")
message(" ")

# Used in report_gen_p2.R
if(!require(seqinr)) {
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}
library(seqinr)

# Used in report_gen_p2.R
if(!require(VariantAnnotation)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("VariantAnnotation")
}
library(VariantAnnotation)

# Used in getDepthStats-parallel.R
if(!require(doParallel)) {
  install.packages('doParallel', repos='http://cran.us.r-project.org')
}  
library(doParallel)  

message("Finished installing R packages used in the pipeline.")
message("****************************************************")
message(" ")

