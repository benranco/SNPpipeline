args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
file <- args[2]
haploidOrDiploid <- args[3]

#path <- "/home/gosuzombie/Desktop/region_38/run"
#file <- "report_p432.Rds"

options(stringsAsFactors = FALSE, warn = 1)

report <- NULL
s <- NULL


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
processTViewOut <- function(tview, rnum, cnum, ref){

    isIndel <- (nchar(ref) > 1)
    
    # get indexes of all * characters in the ref line
    starsInRefLine <- gregexpr(fixed=TRUE, pattern ='*',tview[2,1]) 
    starsInRefLine <- starsInRefLine[[1]][1:length(starsInRefLine[[1]])] # now it's a simple vector
    
    # if there are actually stars in the ref line, update the ref to include those stars:
    if(starsInRefLine[1] != -1) {
      # for each * that is within the ref scope, insert it into the actual ref, so the ref will be the proper
      # length for comparing against the line below it. 
      # (this loop assumes the gregexpr function call above returned the indexes of star characters in 
      # increasing order)  
      for (n in 1:length(starsInRefLine)) {
        starIndex <- starsInRefLine[n]
        if (starIndex <= nchar(ref)) {
          ref <- paste0(substr(ref, 1, starIndex-1), "*", substr(ref, starIndex, nchar(ref)))
        } else {
          break
        }
      }
    }
    
    if(substring(tview[2,1], 1 ,nchar(ref)) == ref)
    {
      if(nrow(tview) >= 3)
      {
        # The following logic is based on information from these two links: 
        # https://en.wikipedia.org/wiki/Pileup_format
        # https://en.wikipedia.org/wiki/Nucleic_acid_sequence
        # It assumes the very first line returned from tview is either some sort of header or blank spaces, 
        # the second line begins with the REF, and the third line is a summary/conclusion of all the lines 
        # below it.
        #
        # pseudocode:
        # only proceed if REF matches the first part of line 2 (first data line)
        # now, for comparing line 3 to line 2:
        # for each position:
        # if . or , then use REF/REF (diploid) or REF/ (haploid)
        # if (A|G|C|T|a|g|c|t) then use uppercase(x/x or x/)
        # if diploid:
        #   if W then A/T or T/A (always REF/ALT, or if neither are REF then the order doesn't matter). 
        #     (haploid shouldn't have these) (if it is part of an indel then put REF before the / and ALT 
        #     after the /) (ACTUALLY, if indel, we don't treat the indels that have WSMKRY, too complex)
        #   if S then C/G or G/C
        #   if M then A/C or C/A
        #   if K then G/T or T/G
        #   if R then A/G or G/A
        #   if Y then C/T or T/C
        # if empty (as is usually the case when the REF position above it is *), or *, skip this position.
        #     (in the case where it is a character but the REF position above it is a *, we just treat it 
        #     with the same rules as anything else, so the result could end up being longer than the REF (but
        #     we've already adjusted the REF that we're comparing it against to be the same length)).
        #     eg. (with REF = TAAAC): 
        #     T*AAAC*AAAAGGGAGATTTTGGA*T*A**G*G*G*G*G**T*A**GTAA***CC*AACAA*T*G*C*CC*C*C*GCACG
        #     . .... ..........C.Y.... Y .  . . . . .  Y .  ....   .. .R... Y K . .. . . S.R..
        #     This case would be: TAAAC/TAAAC
        #
        #     T*AAAC*AAAAGGGAGATTTTGGA*T*A**G*G*G*G*G**T*A**GTAA***CC*AACAA*T*G*C*CC*C*C*GCACG
        #     .A.... ..........C.Y.... Y .  . . . . .  Y .  ....   .. .R... Y K . .. . . S.R..
        #     This case would be: TAAAAC/TAAAAC
        
        firstHalf <- ""
        secondHalf <- ""
        
        abort <- FALSE
        
        for (n in 1:nchar(ref)) {
          x <- substring(tview[3,1], n ,n) # the character we're looking at in this iteration
          
          # NOTE: there is no code to deal with the case when a character in line 3 is " " (or "*") because
          # we just skip it if so.
          
          # 1 == haploid, 2 == diploid
          if (haploidOrDiploid == 1) { # haploid
            if(x == "." || x == ",") {
              firstHalf <- paste0(firstHalf, substring(ref, n ,n) ) 
            } else if(x == "A" || x == "G" || x == "C" || x == "T") {
              firstHalf <- paste0(firstHalf, x) 
            } else if(x == "a" || x == "g" || x == "c" || x == "t") {
              firstHalf <- paste0(firstHalf, toupper(x)) 
            } else if(x == "W" || x == "S" || x == "M" || x == "K" || x == "R" || x == "Y") {
              abort <- TRUE
              break;
            } else if(x != " ") {
              abort <- TRUE
              break;
            }
          } else { # diploid
            if(x == "." || x == ",") {
              r <- substring(ref, n ,n)
              firstHalf <- paste0(firstHalf, r) 
              secondHalf <- paste0(secondHalf, r) 
            } else if(x == "A" || x == "G" || x == "C" || x == "T") {
              firstHalf <- paste0(firstHalf, x) 
              secondHalf <- paste0(secondHalf, x) 
            } else if(x == "a" || x == "g" || x == "c" || x == "t") {
              firstHalf <- paste0(firstHalf, toupper(x)) 
              secondHalf <- paste0(secondHalf, toupper(x)) 
            } else if(isIndel && (x == "W" || x == "S" || x == "M" || x == "K" || x == "R" || x == "Y") ) {
              abort <- TRUE  # we are leaving indels containing W,S,M,K,R,Y as NA.
              break;  
            } else if(x == "W") {
              if (substring(ref, n ,n) == "T") {
                firstHalf <- paste0(firstHalf, "T") 
                secondHalf <- paste0(secondHalf, "A") 
              } else { 
                firstHalf <- paste0(firstHalf, "A") 
                secondHalf <- paste0(secondHalf, "T") 
              }
            } else if(x == "S") {
              if (substring(ref, n ,n) == "G") {
                firstHalf <- paste0(firstHalf, "G") 
                secondHalf <- paste0(secondHalf, "C") 
              } else { 
                firstHalf <- paste0(firstHalf, "C") 
                secondHalf <- paste0(secondHalf, "G") 
              }
            } else if(x == "M") {
              if (substring(ref, n ,n) == "C") {
                firstHalf <- paste0(firstHalf, "C") 
                secondHalf <- paste0(secondHalf, "A") 
              } else { 
                firstHalf <- paste0(firstHalf, "A") 
                secondHalf <- paste0(secondHalf, "C") 
              }
            } else if(x == "K") {
              if (substring(ref, n ,n) == "T") {
                firstHalf <- paste0(firstHalf, "T") 
                secondHalf <- paste0(secondHalf, "G") 
              } else { 
                firstHalf <- paste0(firstHalf, "G") 
                secondHalf <- paste0(secondHalf, "T") 
              }
            } else if(x == "R") {
              if (substring(ref, n ,n) == "G") {
                firstHalf <- paste0(firstHalf, "G") 
                secondHalf <- paste0(secondHalf, "A") 
              } else { 
                firstHalf <- paste0(firstHalf, "A") 
                secondHalf <- paste0(secondHalf, "G") 
              }
            } else if(x == "Y") {
              if (substring(ref, n ,n) == "T") {
                firstHalf <- paste0(firstHalf, "T") 
                secondHalf <- paste0(secondHalf, "C") 
              } else { 
                firstHalf <- paste0(firstHalf, "C") 
                secondHalf <- paste0(secondHalf, "T") 
              }
            } else if(x != " ") {
              abort <- TRUE
              break;
            }
          } 

        } # end for-loop


        if (!abort) {
          # 1 == haploid, 2 == diploid. If it's haploid, we follow the format in the .tab file of "A/",
          # whereas if it's diploid we follow the format in the .tab file of "A/A".
          
          # Update 2022-04-04: Fixed a bug by using <<-, previously was using <- which doesn't assign to a global variable.
          if (haploidOrDiploid == 1) {
            report[rnum,cnum] <<- paste0(firstHalf, "/")
          } else {
            report[rnum,cnum] <<- paste0(firstHalf, "/", secondHalf)
          }
        }

      }
    }

} # end function


#################################################

testProcessTViewOut <- function() {

    
    # Create this test report: 
    #       CHROM POS    REF sample1 sample2
    #          r1  10   AATG      NA      NA
    #          r2  20      A      NA      NA
    #          r3  30      G      NA      NA
    #          r4  40 AATATC      NA      NA
    myReport <- data.frame(c("r1","r2","r3","r4"), c(10,20,30,40), c("AATG","A","G","AATATC"), character(4), character(4))
    names(myReport) <- c("CHROM","POS","REF","sample1","sample2")
    myReport[,4:5] <- NA
    
    # save state of global variables before appropriating them:    
    realReport <- report
    report <<- myReport    
    realHaploidOrDiploid <- haploidOrDiploid
    haploidOrDiploid <<- 2    
    realS <- s
    s <<- 4
        
    # Test indel with Y, this should be NA:
    data1 <- as.matrix(c("151       161       171       181       191       201       211       221       ",
              "AATGAATTTCCACATGCCTTTGAATCTACTTCTATGCTCACTTATGGCATTGGGAGTTTGGACGGGTGTTGGGAAGGAGA",
              "G.Y..............A....C....R......YK............G........................R......",
              "G.....*..........A....C....G......C.............G..............................."))
    processTViewOut(data1, 1, 4, "AATG")
    write(paste0("This should be TRUE: ", is.na(report[1,4])), stdout())
    
    
    # If REF is just A, this should be A/G, if REF is just G, this should be G/A:
    data2 <- as.matrix(c("151       161       171       181       191       201       211       221       ",
              "AATGAATTTCCACATGCCTTTGAATCTACTTCTATGCTCACTTATGGCATTGGGAGTTTGGACGGGTGTTGGGAAGGAGA",
              "R.Y..............A....C....R......YK............G........................R......",
              "G.....*..........A....C....G......C.............G..............................."))
    processTViewOut(data2, 2, 4, "A")
    write(paste0("This should be TRUE: ", "A/G" == report[2,4]), stdout())
    processTViewOut(data2, 3, 4, "G")
    # false because it doesn't proceed if the REF doesn't match the ref in the data line
    write(paste0("This should be TRUE: ", is.na(report[3,4])), stdout()) 
    
    # If REF is just A, this should be A/G, if REF is just G, this should be G/A:
    data3 <- as.matrix(c("151       161       171       181       191       201       211       221       ",
              "GATGAATTTCCACATGCCTTTGAATCTACTTCTATGCTCACTTATGGCATTGGGAGTTTGGACGGGTGTTGGGAAGGAGA",
              "R.Y..............A....C....R......YK............G........................R......",
              "G.....*..........A....C....G......C.............G..............................."))
    processTViewOut(data3, 2, 5, "G")
    write(paste0("This should be TRUE: ", "G/A" == report[2,5]), stdout())
    
    # Test complex indel, If REF is AATATC, this ALT should be GACAATC/GACAATC:
    data4 <- as.matrix(c("151       161       171       181       191       201       211       221       ",
              "AAT*A***TCCACATGCCTTTGAATCTACTTCTATGCTCACTTATGGCATTGGGAGTTTGGACGGGTGTTGGGAAGGAGA",
              "G.C .A  .........A....C....R......YK............G........................R......",
              "G.....*..........A....C....G......C.............G..............................."))
    processTViewOut(data4, 4, 4, "AATATC")
    write(paste0("This should be TRUE: ", "GACAATC/GACAATC" == report[4,4]), stdout())
    
    # Testing for an unsupported character "U", this should give NA:
    data5 <- as.matrix(c("151       161       171       181       191       201       211       221       ",
              "AATGAATTTCCACATGCCTTTGAATCTACTTCTATGCTCACTTATGGCATTGGGAGTTTGGACGGGTGTTGGGAAGGAGA",
              "U.Y..............A....C....R......YK............G........................R......",
              "G.....*..........A....C....G......C.............G..............................."))
    processTViewOut(data5, 3, 5, "A")
    write(paste0("This should be TRUE: ", is.na(report[3,5])), stdout())
    
    #report
    
    message("-----------------------")
    # Now doing the same tests but for haploid:
    report[,4:5] <<- NA
    haploidOrDiploid <<- 1 

    processTViewOut(data1, 1, 4, "AATG")
    write(paste0("This should be TRUE: ", is.na(report[1,4])), stdout())
    
    # extra test:
    processTViewOut(data1, 1, 5, "AA")
    write(paste0("This should be TRUE: ", "GA/" == report[1,5]), stdout())
    
    processTViewOut(data2, 2, 4, "A")
    write(paste0("This should be TRUE: ", is.na(report[2,4])), stdout())
    processTViewOut(data2, 3, 4, "G")
    # false because it doesn't proceed if the REF doesn't match the ref in the data line
    write(paste0("This should be TRUE: ", is.na(report[3,4])), stdout()) 
    
    processTViewOut(data3, 2, 5, "G")
    write(paste0("This should be TRUE: ", is.na(report[2,5])), stdout())
    
    processTViewOut(data4, 4, 4, "AATATC")
    write(paste0("This should be TRUE: ", "GACAATC/" == report[4,4]), stdout())
    
    # extra test:
    processTViewOut(data4, 4, 5, "A")
    write(paste0("This should be TRUE: ", "G/" == report[4,5]), stdout())

    processTViewOut(data5, 3, 5, "A")
    write(paste0("This should be TRUE: ", is.na(report[3,5])), stdout())
        
    #report
    
    # set global variables back to whatever real data may have been in it:
    report <<- realReport
    haploidOrDiploid <<- realHaploidOrDiploid
    s <<- realS
}


#################################################
# Execution begins here (above are functions):

# To test the processTViewOut function, just run source('parallel_process.R') from within an
# R session, it'll give you an error from the below code, but ignore it, then run testProcessTViewOut().

message(paste0("parallel_process: looking up NA data on ", file))
report <- readRDS(paste0(path, "/reporttemp/", file))

if(!("COMBINED" %in% colnames(report)))
{
  s <- 4
}else 
{
  s <- 5
}



for(a in 1:nrow(report))
{
  reference <- report[a, "REF"]


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
        
        tryCatch(
        {              
          # Update 2022-04-04: added header=FALSE to read.delim to deal with the R error "duplicate 'row.names'
          # are not allowed" in the cases where there is a merely whitespace first line instead of a proper
          # header". 
          # This means that the first line will always be the header (even if the header was just whitespace), and
          # the second line will always be the one containing the REF, and the third line will always be the 
          # summary line of the lines below it.
          tviewOut <- as.matrix(read.delim(pipe(cmd), header=FALSE, sep = "\n"))

          processTViewOut(tviewOut, a, b, reference)
          
          #return(NA) # only necessary when doing tryCatch in a function
        },
        error=function(error_message) {
          printErr(cmd, error_message)
          #return(NA) # only necessary when doing tryCatch in a function
        }
        ) # end tryCatch

        
      } 
    } # end inner for-loop
  }
} # end outer for-loop

saveRDS(report, file = paste0(path, "/reporttemp/", substr(file, 0, nchar(file) - 4), "_filled.Rds"))
