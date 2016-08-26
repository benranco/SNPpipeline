args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
options(stringsAsFactors = FALSE, warn = 1)

message("generating subset of report")
report <- read.csv(paste0(path, "/reports/report.csv"))
report <- report[, -1]

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