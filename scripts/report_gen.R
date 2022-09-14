args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
report <- data.frame()
single <- FALSE

message("generating initial report")

#path <- "~/Desktop/SNPpipeline"
options(stringsAsFactors = FALSE, warn = 1)

message(path)
for(pooled in list.files(paste0(path, "/output", "/pooled"))) {
  #message(pooled)
  if(file.info(paste0(path, "/output", "/pooled", "/", pooled))$size > 0) {
    report <- read.delim(paste0(path, "/output", "/pooled", "/", pooled))
  } else {
    report <- data.frame()
  }
}

message("read pooled data")

if(nrow(report) != 0 && ncol(report) != 0) {
  colnames(report)[1] <- "CHROM"
  colnames(report)[4] <- "COMBINED"
} else {
  single <- TRUE
}

message("reading single data and compiling")
for(single1 in list.files(paste0(path, "/output", "/single"))) {
  if(file.info(paste0(path, "/output", "/single", "/", single1))$size > 0) {
    sr <- read.delim(paste0(path, "/output", "/single", "/", single1))
    colnames(sr)[1] <- "CHROM"
    colnames(sr)[4] <- single1
    if(nrow(report) == 0 && ncol(report) == 0) {
      report <- sr
    } else {
      if(single == TRUE) {
        report <- merge(report, sr, by = c("CHROM", "POS", "REF"), all = TRUE)
      } else {
        report <- merge(report, sr, by = c("CHROM", "POS", "REF"), all.x = TRUE)
      }
    }
    #message(single1)
  }
}

write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "report.csv", sep = "/"), row.names=FALSE)

message("report_gen completed")
