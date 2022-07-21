args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
options(stringsAsFactors = FALSE, warn = 1)

message("generating subset of report")
report <- read.csv(paste(path.expand(path), "reports/report.csv", sep = "/"), check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

counter <- 1
filecounter <- 1
split <- 999

numSubsets <- trunc(nrow(report)/split)+1
numDigits <- floor(log10(numSubsets)) + 1 # this just counts the number of digits in a positive integer

while(counter <= nrow(report)) {
  spl <- data.frame()
  if(counter + split > nrow(report)) {
    spl <- report[c(counter:nrow(report)), ]
  }
  else  {
    spl <- report[c(counter:(counter + split)), ]
  }
  saveRDS(spl, file = paste0(path.expand(path), "/reporttemp/report_p",  formatC(filecounter,width=numDigits,flag=0), ".Rds"))
  
  counter <- counter + split + 1
  filecounter <- filecounter + 1
}

message("subset generation complete")
