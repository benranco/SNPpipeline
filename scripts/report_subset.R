args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
options(stringsAsFactors = FALSE, warn = 1)

message("generating subset of report")
report <- read.csv(paste(path.expand(path), "reports/report.csv", sep = "/"), check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

counter <- 1
filecounter <- 1
split <- 999

while(counter <= nrow(report))
{
  spl <- data.frame()
  if(counter + split > nrow(report))
  {
    spl <- report[c(counter:nrow(report)), ]
  }
  else
  {
    spl <- report[c(counter:(counter + split)), ]
  }
  saveRDS(spl, file = paste0(path.expand(path), "/reporttemp/report_p", filecounter, ".Rds"))
  
  counter <- counter + split + 1
  filecounter <- filecounter + 1
}

message("subset generation complete")
