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

# #######################################################################
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

# #######################################################################
message("editing for cutoffs")
  
report <- as.data.frame(report)

if((ncol(report) > 24 && "COMBINED" %in% colnames(report)) || (ncol(report) > 23 && !("COMBINED" %in% colnames(report))))
{
  doFullEditingForCutoffs <- TRUE
}else
{
  doFullEditingForCutoffs <- FALSE
  message("not enough samples for full cutoff editing, just removing non-biallelic rows")
}


isRowNonBiallelic <- function(namesOfGenotypesInRow)
{
  alleleIndex <- 1
  individualAlleles <- character(length(namesOfGenotypesInRow) * 2) # because, for example "A/A" has two alleles, "A/G" has two...
  for ( type in namesOfGenotypesInRow )
  {
    charsInType <- strsplit(type, "")[[1]]
    individualAlleles[alleleIndex] <- charsInType[1]
    alleleIndex <- alleleIndex + 1
    individualAlleles[alleleIndex] <- ifelse (length(charsInType) > 2, charsInType[3], charsInType[1]) # There used to be e.g. "A/" as a shorthand for "A/A"
    alleleIndex <- alleleIndex + 1
  }
  # if there are 3 or more unique alleles in the row, return TRUE, else return FALSE
  return ( ifelse( length(names(table(individualAlleles, useNA = "no"))) >= 3, TRUE, FALSE ) )
}


numrows <- nrow(report)
numColsInReport <- ncol(report)
rowsToRemove <- NULL
rowsToRemove <- numeric(numrows)


if(!("COMBINED" %in% colnames(report)))
{
  startCol <- 4
}else
{
  startCol <- 5
}

numDataCols <- numColsInReport - (startCol - 1)

reportIndex <- 1
rowsToRemoveIndex <- 1

while(reportIndex <= numrows)
{
  rowData <- as.matrix(report[reportIndex,startCol:numColsInReport]) 
  tabulatedRowData <- table(rowData, useNA = "no")
  tabulatedRowDataNames <- names(tabulatedRowData)

  # remove row if there's only one (or 0) type of non-NA data
  if( length(tabulatedRowDataNames) <= 1 & doFullEditingForCutoffs)
  {
    rowsToRemove[rowsToRemoveIndex] <- reportIndex
    rowsToRemoveIndex <- rowsToRemoveIndex + 1
  }
  # remove row if more than half of it's data values are NA
  else if( sum(tabulatedRowData) < numDataCols %/% 2 & doFullEditingForCutoffs)
  {
    rowsToRemove[rowsToRemoveIndex] <- reportIndex
    rowsToRemoveIndex <- rowsToRemoveIndex + 1
  }
  # remove row if it has 3 or more unique alleles (eg. a row with A/A,A/C,C/C is kept, but 
  # a row with A/A,A/C,C/T is removed). We only keep biallelic data.
  else if( isRowNonBiallelic(tabulatedRowDataNames) )
  {
    rowsToRemove[rowsToRemoveIndex] <- reportIndex
    rowsToRemoveIndex <- rowsToRemoveIndex + 1
  }

  reportIndex <- reportIndex + 1
}

rowsToRemove <- rowsToRemove[!(rowsToRemove %in% c(0))] # remove all elements containing 0 from rowsToRemove
if(length(rowsToRemove) > 0) 
{
  report <- report[-rowsToRemove[1:length(rowsToRemove)], ]
}
write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "edited_report.csv", sep = "/"))


# #######################################################################
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

numColsInReport <- ncol(report)
snpp <- as.data.frame(snpp)
numRows <- nrow(snpp)

snpp$A <- integer(numRows)
snpp$C <- integer(numRows)
snpp$T <- integer(numRows)
snpp$G <- integer(numRows)
snpp$empty <- integer(numRows)
snpp$max <- integer(numRows)
snpp$second_max <- integer(numRows)
snpp$sum <- integer(numRows)
snpp$MAF <- numeric(numRows)

for(a in 1:nrow(report))
{
  rowAsMatrix <- as.matrix(report[a, s:numColsInReport])
  factoredGenotypes <- table(rowAsMatrix, useNA = "no")
  totalAlleleCount <- sum(factoredGenotypes) * 2 # because, for example "A/A" has two alleles, "A/G" has two...
  
  # count how many times each allele occurs (factoredGenotypes already lists how many times each genotype occurs)
  for ( type in names(factoredGenotypes) )
  {
    charsInType <- strsplit(type, "")[[1]]
    firstChar <- charsInType[1]
    thirdChar <- ifelse (length(charsInType) > 2, charsInType[3], charsInType[1]) # There used to be e.g. "A/" as a shorthand for "A/A", so count it twice if it still happens to be that way
    snpp[a,firstChar] <- snpp[a,firstChar][[1]] + factoredGenotypes[type][[1]]
    snpp[a,thirdChar] <- snpp[a,thirdChar][[1]] + factoredGenotypes[type][[1]]
  }

  snpp[a, "empty"] <- sum(is.na(rowAsMatrix))
  snpp[a, "max"] <- max(snpp[a, c("A", "C", "T", "G")])
  snpp[a, "second_max"] <- sort(snpp[a, c("A", "C", "T", "G")], TRUE)[2]
  snpp[a, "sum"] <- sum(snpp[a, c("A", "C", "T", "G")])
  snpp[a, "MAF"] <- snpp[a, "second_max"] / (snpp[a, "second_max"] + snpp[a, "max"]) 
  #snpp[a, "chi"] <- ((snpp[a,"max"] - snpp[a,"sum"]%/%2)^2)/snpp[a, "sum"]%/%2  + ((snpp[a,"second_max"] - snpp[a,"sum"]%/%2)^2)/snpp[a, "sum"]%/%2 
}

write.csv(snpp, paste(paste(path.expand(path), "reports", sep = "/"), "percentage_snps.csv", sep = "/"))

# #######################################################################
message("generating MAF_cutoff_report")

snpp <- snpp[order(-snpp$MAF), ]
snpp <- snpp[snpp$MAF >= MAF_CUTOFF,]
report <- report[rownames(snpp), ]

write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report.csv", sep = "/"))

# #######################################################################
message("generating site mutation percentage data")

fastaRef <- read.fasta(file = paste0(path, "/reference/formatted_output.fasta"), as.string = TRUE)
mutationReport <- data.frame()

for(sector in 1:length(names(fastaRef)))
{
  sectorName <- attributes(fastaRef[[sector]])$name
  numRowsForSectorName <- length(grep(sectorName, report[, "CHROM"], fixed=TRUE))
  numCharsInSequence <- nchar(fastaRef[[sector]][1])

  mutationReport[sectorName, "role"] <- attributes(fastaRef[[sector]])$Annot  
  mutationReport[sectorName, "snp"] <- numRowsForSectorName
  mutationReport[sectorName, "length"] <- numCharsInSequence
  mutationReport[sectorName, "percentage SNP"] <- numRowsForSectorName / numCharsInSequence
}

mutationReport <- mutationReport[order(-mutationReport$`percentage SNP`), ]
mutationReport$`percentage SNP` <- mutationReport$`percentage SNP` * 100

write.csv(mutationReport, paste(paste(path.expand(path), "reports", sep = "/"), "mutation_percentage.csv", sep = "/"))

# #######################################################################
message("replacing alleles with characters for chi square test")

# New logic:
# remove rows with NA in MORE than 20% of sites
# calculate occurences of each character (eg. A/A, A/G, G/G = 3 A, 3 G)
#   if the most frequent character has more than 95% of occurences, remove row
#   if the second most frequent character has less than 5% of occurences, remove row (TODO: Or do we sum all the minor characters to get 5%, if there's more than one minor character? - I doubt this, because H,A,B encoding only allows for three characters.)
#   set the heterozygous genotype (eg. G/T) to H, the most frequent homo type to A (eg. G/G), and the second most frequent (if it exists, eg. T/T) to B (if there's a tie in homo type frequency then set the REF/REF one to A)
#   if there are no heterozygous genotypes, set the most frequent homo type to H and the second most to A (if there's a tie in homo type frequency then set the REF/REF one to H)
#   set all samples that contain at least one other allele besides the top two to NA (TODO: is this right?)

# TODO: The above logic is fine for biallelic sites, but what about triallelic sites? They would have three
# possible heterozygous genotypes (REF/ALT1, REF/ALT2, ALT1/ALT2); if I set the most frequently occuring one to
# H, what do I set the other two to? They would also have three possible homozygous genotypes (REF/REF, ALT1/ALT1, ALT2/ALT2); if I set the most frequently occuring one to A and the second most to B, what do I set the other to?

reportc <- report
curRow <- 1
totalRows <- nrow(reportc)

numDataCols <- ncol(reportc) - (s-1)

while(curRow <= totalRows)
{
  datap <- reportc[curRow, s:ncol(reportc)]

  numNAsInRow <- sum(is.na(datap))
  if( numNAsInRow/numDataCols > 0.2 ) 
  {
    # Remove the row if it has more than 20% NA's
    reportc <- reportc[-curRow, ]
    totalRows <- totalRows - 1
    curRow <- curRow - 1
  }
  else
  {
    va <- sort(table(as.matrix(datap), useNA = "no"), decreasing = TRUE)
    totalAlleleCount <- sum(va) * 2 # because, for example "A/A" has two alleles, "A/G" has two...
    alleleSums <- data.frame(A = 0, C = 0, G = 0, T = 0)
    
    # count how many times each allele occurs    
    for ( type in names(va) )
    {
      charsInType <- strsplit(type, "")[[1]]
      firstChar <- charsInType[1]
      thirdChar <- ifelse (length(charsInType) > 2, charsInType[3], charsInType[1]) # There used to be e.g. "A/" as a shorthand for "A/A"
      alleleSums[firstChar] <- alleleSums[firstChar][[1]] + va[type][[1]]
      alleleSums[thirdChar] <- alleleSums[thirdChar][[1]] + va[type][[1]]
    }
    alleleSums <- sort(alleleSums, decreasing = TRUE)
    
    if ( alleleSums[1]/totalAlleleCount > 0.95 )
    {
      # Remove the row if the most frequent allele occurs more than 95% of the time
      reportc <- reportc[-curRow, ]
      totalRows <- totalRows - 1
      curRow <- curRow - 1
    }
    else if ( alleleSums[2]/totalAlleleCount < 0.05 )
    {
      # Remove the row if the second most frequent allele occurs less than 5% of the time
      reportc <- reportc[-curRow, ]
      totalRows <- totalRows - 1
      curRow <- curRow - 1
    }
    
    HType <- NA
    AType <- NA
    BType <- NA

    # figure out which way around the heterozygous genotype should be (eg. "G/T" or "T/G") 
    # by taking the one with the most occurences (va is sorted by number of occurences)
    # (most likely only one of them will have any occurences at all)
    HTypeOpt1 <- paste0(names(alleleSums)[1],"/",names(alleleSums)[2])
    HTypeOpt2 <- paste0(names(alleleSums)[2],"/",names(alleleSums)[1])
    HTypeOpt1Frequency <- 0
    HTypeOpt2Frequency <- 0
    for (type in names(va))
    {
      if (type == HTypeOpt1)
      {
        if (is.na(HType)) 
        { 
          HType <- HTypeOpt1 
        }
        HTypeOpt1Frequency <- va[type]
      }
      else if (type == HTypeOpt2)
      {
        if (is.na(HType)) 
        { 
          HType <- HTypeOpt2
        }
        HTypeOpt2Frequency <- va[type]
      }
    }
    # if HTypeOpt1Frequency and HTypeOpt2Frequency are tied, pick the one that starts with the REF
    if ( HTypeOpt1Frequency > 0 & HTypeOpt1Frequency == HTypeOpt2Frequency )
    {
      if ( names(alleleSums)[1] == reportc[curRow, "REF"] )
      {
        HType <- HTypeOpt1
      }
      else if ( names(alleleSums)[2] == reportc[curRow, "REF"] )
      {
        HType <- HTypeOpt2
      }
    }
    
    firstAlleleIndex <- 1
    secondAlleleIndex <- 2
    # if the first two alleles are tied in terms of frequency, choose the one that matches the REF
    # as the first one.    
    if ( alleleSums[1] == alleleSums[2] & names(alleleSums)[2] == reportc[curRow, "REF"] ) 
    {
      firstAlleleIndex <- 2
      secondAlleleIndex <- 1
    }
    
    # if there was no heterozygous site containing the two most frequent alleles
    if (is.na(HType))
    {
      HType <- paste0(names(alleleSums)[firstAlleleIndex],"/",names(alleleSums)[firstAlleleIndex])
      AType <- paste0(names(alleleSums)[secondAlleleIndex],"/",names(alleleSums)[secondAlleleIndex])
    }
    else
    {
      AType <- paste0(names(alleleSums)[firstAlleleIndex],"/",names(alleleSums)[firstAlleleIndex])
      BType <- paste0(names(alleleSums)[secondAlleleIndex],"/",names(alleleSums)[secondAlleleIndex])
    }
    
    reportc[curRow, s:ncol(reportc)] <- gsub(HType, "H", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = TRUE)
    reportc[curRow, s:ncol(reportc)] <- gsub(AType, "A", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = TRUE)
    if (!is.na(BType))
    {      
      reportc[curRow, s:ncol(reportc)] <- gsub(BType, "B", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = TRUE)
    }
    # replace all elements consisting of a character followed by a "/" followed by an optional character(s) with NA (i.e. all elements that haven't already been replaced with H, A or B):
    reportc[curRow, s:ncol(reportc)] <- gsub("^./.*$", NA, as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = FALSE) 
  }

  curRow <- curRow + 1
}

write.csv(reportc, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report_chi.csv", sep = "/"))

# #######################################################################
message("converting the chi square report to .linkage format")

# Information on .linkage file format (for use with Le-MAP2 software):
#
# Data is tab-delimited.
#
# First 6 Columns contain "pedigree information":
# 1. Family ID (can be alphanumeric)
# 2. Individual ID (must be unique within family, can be alphanumeric)
# 3. Father ID (0 if father is not in family)
# 4. Mother ID (0 if mother is not in family)
# 5. Gender (0=unknown, 1=male, 2=female)
# 6. Affection status (0=unknown, 1=unaffected, 2=affected)
#
# Columns 7 and onward describe the phenotype data, separated by tabs. There 
# are four different types of phenotype data supported by the LINKAGE format 
# (Numbered Alleles, Binary Factors, Affection Status, Quantitative Traits), 
# and the Lep-MAP2 documentation uses Numbered Alleles.
#
# With Numbered Alleles, each genotype is represented as a pair of numbers 
# (eg: 1 2     2 2   1 1    1 2). Each number represents an allele, and 0 
# represents an unknown allele.
#
# For our purposes (input for the Lep-MAP2 software), we are setting:
# H = "1 2", A = "1 1", B = "2 2", NA = "0 0".
#
# Lep-MAP2 documentation: https://sourceforge.net/p/lepmap2/wiki/browse_pages/
#
# Official LINKAGE file format documentation available:
# http://www.jurgott.org/linkage/LinkagePC.html#__RefHeading__137_1806185151
# http://www.jurgott.org/linkage/LinkageUser.pdf

reportLinkageGenotypes <- reportc[ , s:ncol(reportc)]
reportLinkageGenotypes[reportLinkageGenotypes=="H"] <- "1 2"
reportLinkageGenotypes[reportLinkageGenotypes=="A"] <- "1 1"
reportLinkageGenotypes[reportLinkageGenotypes=="B"] <- "2 2"
reportLinkageGenotypes[is.na(reportLinkageGenotypes)] <- "0 0"
reportLinkageGenotypes[reportLinkageGenotypes=="-"] <- "0 0" # in case NA "-" has already been substituted with "-"

reportLinkage <- data.frame(family = reportc$CHROM, id = reportc$POS, fatherId = c(0), motherId = c(0), gender = c(0), affectionStatus = c(0), reportLinkageGenotypes)

write.table(reportLinkage, file= paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report_chi.linkage", sep = "/"), append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# #######################################################################
message("generate probability values")

reportd <- report 

startingCol <- ifelse(("COMBINED" %in% colnames(reportd)), s-1, s) # this makes sure we include the COMBINED col in our for loop if it exists

for(x in startingCol:ncol(reportd))
{
  print(paste0("------------col ",x)) 
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
    ind <- NA
    if(!is.na(reportd[y, x]))
    {
      # The only time two characters occur (other than NA) is the shorthand version of "REF/REF": 
      # "REF/" (eg. if REF is A, instead of writing A/A, write A/). In this case, since they both match
      # the REF, set likelihood value to 1. This is actually legacy code, because reports are no longer 
      # generated with the shorthand "REF/" way.
      if(nchar(reportd[y,x]) == 2 & strsplit(reportd[y,x], "")[[1]][1] == reportd[y, "REF"] & strsplit(reportd[y,x], "")[[1]][2] == "/" )
      {
        reportd[y,x] <- 1
      }
      # three characters, eg. "A/A" or "A/T".
      else if(nchar(reportd[y,x]) == 3)  
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
          # if the 1st and 3rd match each other:
          if(strsplit(reportd[y,x], "")[[1]][1] == strsplit(reportd[y,x], "")[[1]][3] )
          {
            # if they match the REF:
            if(strsplit(reportd[y,x], "")[[1]][1] == reportd[y, "REF"] )
            {
              # (GL position 1)
              reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][1]
              # TODO: Not yet sure if I just want to set it to 1 automatically here, but if so, 
              # use this code instead:
              #reportd[y,x] <- 1
            }
            # if they DON'T match the REF:
            else 
            {
              # Check the Genotype (GT) to determine which GL position to use
              if( as.list(ap$`*:*-*`$GENO$GT[ind])[[1]][1] == "1/1" ) # "1/1" means GL position 3
              {
                reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][3]
              }
              else # "2/2" means GL position 6
              {
                reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][6]
              }
            }
          }     
          # if the 1st and 3rd DON'T match each other:
          else
          {
            # if the first character matches the REF:
            if(strsplit(reportd[y,x], "")[[1]][1] == reportd[y, "REF"] )
            {
              # Check the Genotype (GT) todetermine which GL position to use
              if( as.list(ap$`*:*-*`$GENO$GT[ind])[[1]][1] == "0/1" ) # "0/1" means GL position 2
              {
                reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][2]
              }
              else # "0/2" means GL position 4
              {
                reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][4]
              }
            }
            # if the first character DOESN'T match match the REF:
            else 
            {
              # it must be GL position 5 (only possible for triallelic sites)
              reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][5]
            }
          }
        }
        else # if ind is NA or null or has no value
        {
          # if the 1st and 3rd character match each other and match the REF
          if( strsplit(reportd[y,x], "")[[1]][1] == reportd[y, "REF"] & strsplit(reportd[y,x], "")[[1]][1] == strsplit(reportd[y,x], "")[[1]][3])
          {
            reportd[y,x] <- 1
          }
          else
          {
            reportd[y, x] <- NA
          }
        }
      } # end of 3 character cases
      else
      {
        reportd[y, x] <- NA
      }

      # Old Comments:
      # TODO: If we raise 10 ^ (much greater than -320) I think the number is too large, resulting in 0.


    } # end of if(!is.na(reportd[y, x]))
  } # end for loop (iterating down the rows)
} # end for loop (iterating across the colums)

write.csv(reportd, paste(paste(path.expand(path), "reports", sep = "/"), "probability.csv", sep = "/"))


message("report_gen part 2 complete")
