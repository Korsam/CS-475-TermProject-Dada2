#########################################################
###             Example 4: fastq files                ###
###                                                   ###
###      Standard dada2 batch file for myPhyloDB      ###
###       Advanced users, may wish to change the      ###
###   pipeline; however, please note that file names  ###
###       will need to be changed accordingly.        ###
###                                                   ###
###        created for: myPhyloDB v1.2.1              ###
###              www.myphylodb.org                    ###
###        created by: Daniel K. Manter               ###
###              daniel.manter@ars.usda.gov           ###
###        date: August 23, 2017                      ###
#########################################################


startTime <- Sys.time()
################### Check and load for required packages ###############
.cran_packages <-  c("ggplot2", "gridExtra", "Biostrings")
.bioc_packages <- c("dada2", "phyloseq")
sapply(c(.cran_packages, .bioc_packages), require, character.only=TRUE)
checkTime <- Sys.time()
#########################################################################


############################# Get data ##################################
set.seed(100)
path <- ('input')
dir.create('Forward')
dir.create('Reverse')

file_info <- read.csv('temp.files', sep='\t', header=FALSE)
sample.names <- file_info[,1]
sample.names

fastqFs <- file_info[,2]
fastqFs

fastqRs <- file_info[,3]
fastqRs
dataTime <- Sys.time()
#########################################################################



#################### Remove Primers w/ cutadapt ##########################
dir.create('Forward/Cut')
dir.create('Reverse/Cut')
cutpathF <- file.path("Forward/Cut")
cutpathR <- file.path("Reverse/Cut")

# select one forward primer
#primerF <- 'AGAGTTTGATCMTGGCTCAG' # 27F
primerF <- 'CCTACGGGNGGCWGCAG'    # 341f
#primerF <- 'GTGYCAGCMGCCGCGGTAA'  # 515fb

# select one reverse primer
#primerR <- 'TTACCGCGGCKGCTGGCAC' # 515r
primerR <- 'GACTACHVGGGTATCTAATCC'  # 806rb
  
# remove primers
for(i in 1:length(fastqFs)) {
  infileF = file.path(path, fastqFs[i])
  outfileF = file.path(cutpathF, fastqFs[i])
  infileR = file.path(path, fastqRs[i])
  outfileR = file.path(cutpathR, fastqRs[i])
  
  cmd <- paste0('cutadapt',
			   ' -g ^', primerF, ' -G ^', primerR, 
			   ' -o ', outfileF, ' -p ', outfileR, 
			   ' ', infileF, ' ', infileR)
  system(cmd)
}

# create trimmedreference database
infile = 'gg_13_5_99.dkm.fa.gz'
outfile = 'ref_trimmed.fa.gz'

rc_primerR <- as.character(reverseComplement(DNAString(primerR)))
cmd <- paste0('cutadapt --trimmed-only',
              ' -g ', primerF, '...', rc_primerR,
              ' -o ', outfile, 
              ' ', infile)
system(cmd)
cutTime <- Sys.time()
#########################################################################


######################### Trim and filter ###############################
# Filtering: THESE PARAMETERS MAY NOT BE OPTIMAL FOR ALL DATASETS
dir.create('Forward/Filtered')
dir.create('Reverse/Filtered')
filtpathF <- file.path("Forward/Filtered")
filtpathR <- file.path("Reverse/Filtered")

fastqFs <- list.files(cutpathF) # recreate list in case some samples did not have data
fastqRs <- list.files(cutpathR) # recreate list in case some samples did not have data

out <- filterAndTrim(fwd=file.path(cutpathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(cutpathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, 
              multithread=4)
out
trimTime <- Sys.time()
#########################################################################



##################### Infer sequence variants ###########################
filtFs <- list.files(filtpathF, full.names=TRUE)
filtRs <- list.files(filtpathR, full.names=TRUE)

# need to match file list with original sample names due to potential sorting issues
# or missing data
filtFs.base <- as.data.frame(list.files(filtpathF, full.names=FALSE))
names(filtFs.base) <- c('dir_list')
new_file_info <- merge(filtFs.base, file_info, by.x='dir_list', by.y='V2', sort=FALSE)
sample.names <- new_file_info$V1

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Learn error rates
set.seed(100)
# preset work environment gets to here

preErrTime <- Sys.time()
errF <- learnErrors(filtFs, nreads=1e6, multithread=4)
errR <- learnErrors(filtRs, nreads=1e6, multithread=4)
postErrTime <- Sys.time()

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=4)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=4)
  merger <- mergePairs(ddF, derepF, ddR, derepR,
                       minOverlap=20, maxMismatch=0,
                       verbose=TRUE)
  mergers[[sam]] <- merger
}

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab, method="pooled", multithread=4)
seqTime <- Sys.time()
##########################################################################



############################ Export Files ################################
# Representative Sequences
numOtus <- ncol(seqtab)
uniquesToFasta(seqtab, 
  ids=paste0("OTU", seq(numOtus)),
  fout='dada.fasta'  
)

# Mothur shared file
otab <- otu_table(seqtab, taxa_are_rows=FALSE)
colnames(otab) <- paste0("OTU", seq(ncol(otab)))
otab <- data.frame(otab)
df <- cbind(label='dada', group=row.names(otab), numOtus=numOtus, otab)
write.table(df, 
  "dada.shared",
  quote=FALSE, sep="\t", row.names=F
)
expTime <- Sys.time()
##########################################################################

cat("Completed, timestamps to follow")
cat("Start:\t", capture.output(Sys.time()-startTime), "\n")
cat("Check:\t", capture.output(checkTime-startTime), "\n")
cat("Data:\t", capture.output(dataTime-checkTime), "\n")
cat("Cut:\t", capture.output(cutTime-dataTime), "\n")
cat("Trim:\t", capture.output(trimTime-cutTime), "\n")
cat("Seq:\t", capture.output(seqTime-trimTime), "\n")
cat("Export:\t", capture.output(expTime-seqTime), "\n")

cat("learnErrors:\t", capture.output(postErrTime-preErrTime), "\n")
cat("process:\t", capture.output(seqTime-postErrTime), "\n")



