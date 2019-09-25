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

#################################################
printErr <- function(command, err) {
    message("**************** parallel_process.R: An error occured trying to execute the following system command: ")
    message(command)
    message("And below is the error message from R:")
    message(err)
    message(" ")
    message("We will skip trying to retrieve data from the BAM file for this particular row+sample.")
    return(NA)
}

#################################################
searchBAMfile <- function(fn, cmd, rnum, cnum, ref, ref_length){
  tryCatch(
  {        
    out <- as.matrix(read.delim(pipe(cmd), sep = "\n"))
    
    if(substring(out[1,1], 1 ,ref_length) == ref)
    {
      if(nrow(out) >= 2)
      {
        match <- TRUE
        for (n in 1:ref_length) {
          if( !( substring(out[2,1], n ,n) == "." || substring(out[2,1], n ,n) == "," ) )
          {
            match <- FALSE
            break
          }
        } # end for-loop

        if( match == TRUE ) 
        {       
          # 1 == haploid, 2 == diploid. If it's haploid, we follow the format in the .tab file of "A/",
          # whereas if it's diploid we follow the format in the .tab file of "A/A".
          if (haploidOrDiploid == 1) {
            report[rnum,cnum] <- paste0(ref, "/")
          } else {
            report[rnum,cnum] <- paste0(ref, "/", ref)
          }
        }

      }
    }
    return(NA)
  },
  error=function(error_message) {
    printErr(cmd, error_message)
    return(NA)
  }
  ) # end tryCatch
}

#################################################

for(a in 1:nrow(report))
{
  reference <- report[a, "REF"]
  refLength <- nchar(reference)

  if(s <= ncol(report))
  {
    for(b in s:ncol(report))
    {
      if(is.na(report[a,b]))
      {
        fn <- paste(substr(colnames(report)[b], 0 , nchar(colnames(report)[b]) - 4), 
                    "_sorted_markDup.bam", sep = "")
        cmd <- paste0(path, "/tools/samtools-1.3.1/samtools tview ", path ,"/dataTemp/single/", fn ,
                      " ", path, "/reference/formatted_output.fasta -d T", 
                      ' -p \"', paste(report[a, "CHROM"], report[a, "POS"], sep = ":"), '"')
        searchBAMfile(fn, cmd, a, b, reference, refLength)
      } 
    } # end inner for-loop
  }
} # end outer for-loop

saveRDS(report, file = paste0(path, "/reporttemp/", substr(file, 0, nchar(file) - 4), "_filled.Rds"))
