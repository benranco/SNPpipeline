args <- commandArgs(trailingOnly = TRUE)
path <- args[1]

MAF_CUTOFF <- args[2]
MAF_CUTOFF <- as.double(MAF_CUTOFF)

#debug arguments
#path <- "/home/gosuzombie/Desktop/region38"
#MAF_CUTOFF <- 0.3

report <- data.frame()
options(stringsAsFactors = FALSE, warn = 1)

message("running report generation part 2")
if(!require(seqinr))
{
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}

if(!require(VariantAnnotation))
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("VariantAnnotation")
}

library(VariantAnnotation)
library(seqinr)

message("reading and merging split data")
for(single1 in list.files(paste0(path, "/reporttemp")))
{
  if(grepl("filled.Rds", single1))
  {
    #message(single1)
    sr <- readRDS(paste0(path, "/reporttemp/", single1))
    if(nrow(report) == 0 && ncol(report) == 0)
    {
      report <- sr
    }
    else
    {
      report <- rbind(report, sr)
    }
  }
}

write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "filled_report.csv", sep = "/"))

if((ncol(report) > 24 && "COMBINED" %in% colnames(report)) || (ncol(report) > 23 && !("COMBINED" %in% colnames(report))))
{
  message("editing for cutoffs")
  tr <- report
  
  report <- tr
  report <- as.data.frame(report)
  nc <- nrow(report)
  
  a <- 1
  
  if(!("COMBINED" %in% colnames(report)))
  {
    rfc <- 4
  }else
  {
    rfc <- 5
  }
  
  while(a <= nc)
  {
    c <- as.matrix(report[a, c(1:rfc) * -1])
    if(sum(is.na(c)) > length(c)%/%2)
    {
      report <- report[-a, ]
      a <- a -1
      nc <- nc - 1
    }
    else
    {
      if((length(unique(as.factor(c))) == 2 && (NA %in% c)) || (length(unique(as.factor(c))) == 1 && !(NA %in% c))) 
      {
        report <- report[-a, ]
        a <- a -1
        nc <- nc - 1
      }
    }
    a <- a + 1
  }
  write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "edited_report.csv", sep = "/"))
}else
{
  message("not enough samples, skipping cutoff editing")
}

message("finding snp percentage per site")
if(!("COMBINED" %in% colnames(report)))
{
  snpp <- report[, c(1:3)]
  s <- 4
}else
{
  snpp <- report[, c(1:4)]
  s <- 5
}

snpp <- as.data.frame(snpp)
for(a in 1:nrow(report))
{
  snpp[a, "A"] <- length(grep("A", as.matrix(report[a, c(s:length(report[a,]))])))
  snpp[a, "C"] <- length(grep("C", as.matrix(report[a, c(s:length(report[a,]))])))
  snpp[a, "T"] <- length(grep("T", as.matrix(report[a, c(s:length(report[a,]))])))
  snpp[a, "G"] <- length(grep("G", as.matrix(report[a, c(s:length(report[a,]))])))
  snpp[a, "empty"] <- sum(is.na(as.matrix(report[a, c(s:length(report[a,]))])))
  snpp[a, "max"] <- max(snpp[a, c("A", "C", "T", "G")])
  snpp[a, "second_max"] <- sort(snpp[a, c("A", "C", "T", "G")], TRUE)[2]
  snpp[a, "sum"] <- sum(snpp[a, c("A", "C", "T", "G")])
  snpp[a, "MAF"] <- snpp[a, "second_max"] / (snpp[a, "second_max"] + snpp[a, "max"]) 
  #snpp[a, "chi"] <- ((snpp[a,"max"] - snpp[a,"sum"]%/%2)^2)/snpp[a, "sum"]%/%2  + ((snpp[a,"second_max"] - snpp[a,"sum"]%/%2)^2)/snpp[a, "sum"]%/%2 
}

write.csv(snpp, paste(paste(path.expand(path), "reports", sep = "/"), "percentage_snps.csv", sep = "/"))

snpp <- snpp[order(-snpp$MAF), ]
snpp <- snpp[snpp$MAF >= MAF_CUTOFF,]
report <- report[rownames(snpp), ]

message("generating site mutation percentage data")

output <- read.fasta(file = paste0(path, "/reference/formatted_output.fasta"), as.string = TRUE)
out <- data.frame()

for(sector in 1:length(names(output)))
{
  out[attributes(output[[sector]])$name, "role"] <- attributes(output[[sector]])$Annot  
  out[attributes(output[[sector]])$name, "snp"] <- length(grep(attributes(output[[sector]])$name, report[, "CHROM"]))
  out[attributes(output[[sector]])$name, "length"] <- nchar(output[[sector]][1])
  out[attributes(output[[sector]])$name, "percentage SNP"] <- length(grep(attributes(output[[sector]])$name, report[, "CHROM"])) / nchar(output[[sector]][1])
}

out <- out[order(-out$`percentage SNP`), ]
out$`percentage SNP` <- out$`percentage SNP` * 100

write.csv(out, paste(paste(path.expand(path), "reports", sep = "/"), "mutation_percentage.csv", sep = "/"))

write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report.csv", sep = "/"))

message("replacing alleles with characters for chi square test")

reportc <- report
ct <- 1
cl <- nrow(reportc)

while(ct <= cl)
{
  datap <- reportc[ct, s:ncol(reportc)]
  va <- sort(table(as.matrix(datap)), decreasing = TRUE)
  if(length(va) == 1)
  {
    reportc[ct, s:ncol(reportc)] <- gsub(names(va)[1], "H", as.matrix(reportc[ct, s:ncol(reportc)]), fixed = TRUE)
  }
  else if(length(va) > 1 & ((va[1] + va[2])/sum(va) > 0.9))
  {
    # if names(va)[2] contains names(va)[1], do the string substitution for names(va)[2] first so as to not mess up occurences of names(va)[2] by replacing names(va)[1] first
    if(grepl(names(va)[1], names(va)[2]))  
    {
      reportc[ct, s:ncol(reportc)] <- gsub(names(va)[2], "A", as.matrix(reportc[ct, s:ncol(reportc)]), fixed = TRUE)
      reportc[ct, s:ncol(reportc)] <- gsub(names(va)[1], "H", as.matrix(reportc[ct, s:ncol(reportc)]), fixed = TRUE)
    } 
    # otherwise replace names(va)[1] first
    else
    {
      reportc[ct, s:ncol(reportc)] <- gsub(names(va)[1], "H", as.matrix(reportc[ct, s:ncol(reportc)]), fixed = TRUE)
      reportc[ct, s:ncol(reportc)] <- gsub(names(va)[2], "A", as.matrix(reportc[ct, s:ncol(reportc)]), fixed = TRUE)
    }

    if(length(va) > 2)
    {
      for(a in c(3:length(va)))
      {
        reportc[ct, s:ncol(reportc)] <- gsub(names(va)[a], NA , as.matrix(reportc[ct, s:ncol(reportc)]), fixed = TRUE)
      }
    }
  }
  else
  {
    reportc <- reportc[-ct, ]
    cl <- cl - 1
    ct <- ct - 1
  }
  ct <- ct + 1
}

write.csv(reportc, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report_chi.csv", sep = "/"))

message("generate probability values")

reportd <- reportc 

startingCol <- ifelse(("COMBINED" %in% colnames(reportd)), s-1, s)

for(x in startingCol:ncol(reportd))
{
  if(colnames(reportd)[x] == "COMBINED")
  {
    for(cutf in list.files(paste0(path, "/outputTemp/pooled")))
    {
      if(grepl("cutoff", cutf))
      {
        fil <- paste0(path, "/outputTemp/pooled/", cutf)
      }
    }
  }
  else
  {
    fil <- paste0(path, "/outputTemp/single/", substr(colnames(reportd)[x], 0, nchar(colnames(reportd)[x]) - 4), "_cutoff")
  }
  ap <- scanVcf(fil)
  rr <- ap$`*:*-*`$rowRanges
  for(y in 1:nrow(reportd))
  {
    ind <- NA
    if(!is.na(reportd[y, x]))
    {
      if(nchar(reportd[y,x]) == 2)
      {
        if(strsplit(reportd[y,x], "")[[1]][1] != reportd[y, "REF"])
        {
          ind <- intersect(grep(paste0(":", reportd[y, "POS"], "_"), names(rr)), grep(reportd[y, "CHROM"], names(rr)))
          if(!is.null(ind) && length(ind) > 0 && !is.na(ind))
          {
            reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][2]
          }
          else
          {
            reportd[y, x] <- NA
          }
        }
        else
        {
          reportd[y,x] <- 1
        }
      }
      else # added check for >= 3 by Ben (may want to add it back, as well as else clause below)
      {
        ind <- intersect(grep(paste0(":", reportd[y, "POS"], "_"), names(rr)), grep(reportd[y, "CHROM"], names(rr)))        
        if(!is.null(ind) && length(ind) > 0 && !is.na(ind))
        {
          if( (nchar(reportd[y,x]) >= 3) & strsplit(reportd[y,x], "")[[1]][1] != strsplit(reportd[y,x], "")[[1]][3])
          { # added check for >= 3 by Ben
            reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][2]
          }
          else
          {
            reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][3]
          }
        }
        else # added by Ben
        {
          reportd[y, x] <- NA
        }
      }
#      else # added by Ben
#      {
#        reportd[y, x] <- NA
#      }
    }
  }
} 

write.csv(reportd, paste(paste(path.expand(path), "reports", sep = "/"), "probability.csv", sep = "/"))

message("report_gen part 2 complete")
