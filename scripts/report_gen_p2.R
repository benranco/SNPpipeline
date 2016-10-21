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
  if(grepl("filled.Rds", single1, fixed=TRUE))
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
  
  
  
  report <- as.data.frame(report)
  numrows <- nrow(report)
  
  a <- 1
  
  if(!("COMBINED" %in% colnames(report)))
  {
    startCol <- 4
  }else
  {
    startCol <- 5
  }
  
  while(a <= numrows)
  {
    rowData <- as.matrix(report[a,startCol:ncol(report)]) 
    if(sum(is.na(rowData)) > length(rowData)%/%2)
    {
      report <- report[-a, ]
      a <- a -1
      numrows <- numrows - 1
    }
    else
    {
      if((length(unique(as.factor(rowData))) == 2 && (NA %in% rowData)) || (length(unique(as.factor(rowData))) == 1 && !(NA %in% rowData))) 
      {
        report <- report[-a, ]
        a <- a -1
        numrows <- numrows - 1
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
  snpp[a, "A"] <- length(grep("A", as.matrix(report[a, c(s:length(report[a,]))]), fixed=TRUE))
  snpp[a, "C"] <- length(grep("C", as.matrix(report[a, c(s:length(report[a,]))]), fixed=TRUE))
  snpp[a, "T"] <- length(grep("T", as.matrix(report[a, c(s:length(report[a,]))]), fixed=TRUE))
  snpp[a, "G"] <- length(grep("G", as.matrix(report[a, c(s:length(report[a,]))]), fixed=TRUE))
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

fastaRef <- read.fasta(file = paste0(path, "/reference/formatted_output.fasta"), as.string = TRUE)
mutationReport <- data.frame()

for(sector in 1:length(names(fastaRef)))
{
  mutationReport[attributes(fastaRef[[sector]])$name, "role"] <- attributes(fastaRef[[sector]])$Annot  
  mutationReport[attributes(fastaRef[[sector]])$name, "snp"] <- length(grep(attributes(fastaRef[[sector]])$name, report[, "CHROM"], fixed=TRUE))
  mutationReport[attributes(fastaRef[[sector]])$name, "length"] <- nchar(fastaRef[[sector]][1])
  mutationReport[attributes(fastaRef[[sector]])$name, "percentage SNP"] <- length(grep(attributes(fastaRef[[sector]])$name, report[, "CHROM"], fixed=TRUE)) / nchar(fastaRef[[sector]][1])
}

mutationReport <- mutationReport[order(-mutationReport$`percentage SNP`), ]
mutationReport$`percentage SNP` <- mutationReport$`percentage SNP` * 100

write.csv(mutationReport, paste(paste(path.expand(path), "reports", sep = "/"), "mutation_percentage.csv", sep = "/"))

write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report.csv", sep = "/"))

message("replacing alleles with characters for chi square test")

reportc <- report
curRow <- 1
totalRows <- nrow(reportc)

while(curRow <= totalRows)
{
  datap <- reportc[curRow, s:ncol(reportc)]
  va <- sort(table(as.matrix(datap), useNA = "ifany"), decreasing = TRUE) # TODO: Check if useNA is what we want here. Ben

  # TODO: (Fixed, but confirm): What if either of the most frequent values is NA? Do we delete the row? What if there's only two values and one of them is NA? This if statement below will automatically keep the row even if the majority is NA, because NA is not included in va, and therefore is subtracted from the total. I think the table function used above excludes NA's from the factorization/count (I fixed this by adding the useNA parameter in the table function and updating the if statements to include NAs in their calculations, rather than ignoring the NAs and coming up with false calculations because of it).
  if( va[1]/sum(va) > 0.9 & !is.na(names(va)[1]) ) # TODO: Check if !is.na is what we want here. Ben
  {
    reportc[curRow, s:ncol(reportc)] <- gsub(names(va)[1], "H", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = TRUE)
  }
  else if( length(va) > 1 & (va[1]+va[2])/sum(va) > 0.9 & !is.na(names(va)[1]) & !is.na(names(va)[2]) ) # TODO: Check if !is.na is what we want here. Ben
  {
    # if names(va)[2] contains names(va)[1], do the string substitution for names(va)[2] first so as to not mess up occurences of names(va)[2] by replacing names(va)[1] first
    if(grepl(names(va)[1], names(va)[2], fixed=TRUE))  
    {
      reportc[curRow, s:ncol(reportc)] <- gsub(names(va)[2], "A", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = TRUE)
      reportc[curRow, s:ncol(reportc)] <- gsub(names(va)[1], "H", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = TRUE)
    } 
    # otherwise replace names(va)[1] first
    else
    {
      reportc[curRow, s:ncol(reportc)] <- gsub(names(va)[1], "H", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = TRUE)
      reportc[curRow, s:ncol(reportc)] <- gsub(names(va)[2], "A", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = TRUE)
    }

    if(length(va) > 2)
    {
      for(a in c(3:length(va)))
      {
        reportc[curRow, s:ncol(reportc)] <- gsub(names(va)[a], NA , as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = TRUE)
      }
    }
  }
  else
  {
    reportc <- reportc[-curRow, ]
    totalRows <- totalRows - 1
    curRow <- curRow - 1
  }
  curRow <- curRow + 1
}

write.csv(reportc, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report_chi.csv", sep = "/"))

message("generate probability values")

reportd <- report 

startingCol <- ifelse(("COMBINED" %in% colnames(reportd)), s-1, s)

for(x in startingCol:ncol(reportd))
{
#  print("------------") # delete this
  if(colnames(reportd)[x] == "COMBINED")
  {
    for(cutf in list.files(paste0(path, "/outputTemp/pooled")))
    {
      if(grepl("cutoff", cutf, fixed=TRUE))
      {
        fil <- paste0(path, "/outputTemp/pooled/", cutf)
      }
    }
  }
  else
  {
    fil <- paste0(path, "/outputTemp/single/", substr(colnames(reportd)[x], 0, nchar(colnames(reportd)[x]) - 4), "_cutoff")
  }
  # For documentation on the VariantAnnotation packaged used for scanVcf and related operations, see:
  #     http://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html
  # For specification of the vcf file format, see:
  #     https://github.com/samtools/hts-specs (has canonical specifications of VCF format)
  #     http://www.htslib.org/doc/vcf.html (seems to briefly document the VCF format)
  # Of particular relevance to the code below, see the descriptions of GT and GL in GENO, section 1.4.2, 
  # page 5 of:    
  #     http://samtools.github.io/hts-specs/VCFv4.1.pdf
  ap <- scanVcf(fil) # import vcf file (see page 35 of http://bioconductor.org/packages/release/bioc/manuals/VariantAnnotation/man/VariantAnnotation.pdf)
  rr <- ap$`*:*-*`$rowRanges
  for(y in 1:nrow(reportd))
  {
#    tmp = reportd[y,x] # delete this
    ind <- NA
    if(!is.na(reportd[y, x]))
    {
      if(nchar(reportd[y,x]) == 2)
      {
        if(strsplit(reportd[y,x], "")[[1]][1] != reportd[y, "REF"])
        {
          # Set "ind" to the list of the indices of those row range names which contain both this row's "CHROM"
          # value and ":<this row's POS value>_". For our purposes, we only care about the first one in the
          # list if there are more than one. It'll be used to look up the Genotype Likelihood values for the
          # current cell.
          # (Using the fixed=TRUE parameter in grep prevents it from evaluating the first string as a regular
          # expression - important, because it would cause problems for special characters otherwise).
          # (Using the value=FALSE parameter in grep makes it return the integer indices of matching elements,
          # rather than the elements themselves. It's the default, but it doesn't hurt to be specific). 
          ind <- intersect(grep(paste0(":", reportd[y, "POS"], "_"), names(rr), value=FALSE, fixed=TRUE), grep(reportd[y, "CHROM"], names(rr), value=FALSE, fixed=TRUE))
          if(!is.null(ind) && length(ind) > 0 && !is.na(ind))
          {
            # set the cell to 10 ^ (the 2nd of 3 genotype likelihoods in the vcf “_cutoff” file's $GENO$GL[index] record). See  http://samtools.github.io/hts-specs/VCFv4.1.pdf, section 1.4.2, page 5 for info on GL.
            reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][2]
          }
          else
          {
            reportd[y, x] <- NA
          }
        }
        else
        {
#          cat("+") # delete this
          reportd[y,x] <- 1
        }
      }
# TODO: To be uncommented if this is required and correct. (see TODO below). I must first make sure that every time we have a matching sample (eg. G/G) that also matches the REF (eg. G), we can automatically set the cell to 1. Priliminary testing suggests that we can't assume this.
#      else if( (nchar(reportd[y,x]) >= 3) & strsplit(reportd[y,x], "")[[1]][1] == reportd[y, "REF"] & strsplit(reportd[y,x], "")[[1]][1] == strsplit(reportd[y,x], "")[[1]][3])
#      {
##        cat("*") # delete this
#        reportd[y,x] <- 1
#      }
      else if(nchar(reportd[y,x]) >= 3) # added check for >= 3 by Ben 
      {
        # TODO: If I try setting one of the sample values that was just two characters ("G/") to three characters ("G/G"), the ind below evaluates to null or na. It doesn't do this for samples that are naturally three characters with the 1st and 3rd matching. Should I make an extra conditional to check if 1st and 3rd characters match AND they match the REF, and set the new value to 1 if so? This would replicate what is being done for the two character situation above. I've implemented and tested it above, and commented it out for now until I know for sure that this is what it should do. 
        ind <- intersect(grep(paste0(":", reportd[y, "POS"], "_"), names(rr), value=FALSE, fixed=TRUE), grep(reportd[y, "CHROM"], names(rr), value=FALSE, fixed=TRUE))     
        if(!is.null(ind) && length(ind) > 0 && !is.na(ind))
        {
          # TODO: Should we have a case for when the first character matches the REF allele but the third doesn't? Or when the third character matches the REF allele but the first doesn't? And a case for when neither match the REF allele? Would these make any difference?
          if( strsplit(reportd[y,x], "")[[1]][1] != strsplit(reportd[y,x], "")[[1]][3] )
          { # added check for >= 3 by Ben
#            cat(paste0(":",0.0 + as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][3]))  # delete this
            reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][2] # the original setup
            #TODO: Right now, the way these as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][2] are set up, it's always accessing the one that's 0, and raising 10^0 is always 1. I tried switching the 2 and the 3 around. Check to make sure which one is right.
            #TODO: If we raise 10 ^ (much greater than -320) I think the number is too large, resulting in 0.
            # Should we be doing 10 ^ at all?
          }
          else
          {
#            cat(paste0(":",0.0 + as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][2]))  # delete this
            reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][3] # the original setup
          }
        }
        else # added by Ben
        {
#          cat("_") # delete this
          reportd[y, x] <- NA
        }
        
      }
      else # added by Ben
      {
        reportd[y, x] <- NA
      }
    }
#    print(paste0(y, ": REF=", reportd[y, "REF"]," sample: ", tmp, ", result: ",reportd[y,x])) # delete this
  }
} 

write.csv(reportd, paste(paste(path.expand(path), "reports", sep = "/"), "probability.csv", sep = "/"))

message("report_gen part 2 complete")
