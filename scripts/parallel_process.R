args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
file <- args[2]
haploidOrDiploid <- args[3]

#path <- "/home/gosuzombie/Desktop/region_38/run"
#file <- "report_p432.Rds"

#message(paste0("looking up NA data on ", file))
options(stringsAsFactors = FALSE, warn = 1)
message("running parallel process")
report <- readRDS(paste0(path, "/reporttemp/", file))

if(!("COMBINED" %in% colnames(report)))
{
  s <- 4
}else 
{
  s <- 5
}



searchBAMfile <- function(fn, cmd){
  tryCatch(
  {        
    out <- as.matrix(read.delim(pipe(cmd), sep = "\n"))
    
    if(substring(out[1,1], 1 ,1) == report[a, "REF"])
    {
      if(nrow(out) >= 2)
      {
        if(substring(out[2,1], 1 ,1) == "." || substring(out[2,1], 1 ,1) == ",")
        {
          # 1 == haploid, 2 == diploid. If it's haploid, we follow the format in the .tab file of "A/",
          # whereas if it's diploid we follow the format in the .tab file of "A/A".
          if (haploidOrDiploid == 1) {
            report[a,b] <- paste0(report[a, "REF"], "/")
          } else {
            report[a,b] <- paste0(report[a, "REF"], "/", report[a, "REF"])
          }
        }
      }
    }
    return(NA)
  },
  error=function(error_message) {
    message("**************** parallel_process.R: An error occured trying to execute the following system command: ")
    message(cmd)
    message("And below is the error message from R:")
    message(error_message)
    message(" ")
    message("We will skip trying to retrieve data from the BAM file for this particular row+sample.")
    return(NA)
  }
  ) # end tryCatch
}



for(a in 1:nrow(report))
{
  if(s <= ncol(report))
  {
    for(b in s:ncol(report))
    {
      if(is.na(report[a,b]))
      {
        fn <- paste(substr(colnames(report)[b], 0 , nchar(colnames(report)[b]) - 4), 
                    "_sorted.bam", sep = "")
        cmd <- paste0(path, "/tools/samtools-1.3.1/samtools tview ", path ,"/dataTemp/single/", fn ,
                      " ", path, "/reference/formatted_output.fasta -d T", 
                      ' -p \"', paste(report[a, "CHROM"], report[a, "POS"], sep = ":"), '"')
        searchBAMfile(fn, cmd)
      } 
    } # end inner for-loop
  }
} # end outer for-loop

saveRDS(report, file = paste0(path, "/reporttemp/", substr(file, 0, nchar(file) - 4), "_filled.Rds"))
