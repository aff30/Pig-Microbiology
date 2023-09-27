set.seed(200789)

## R script to run DADA2 algorithm on trimmed paired-end reads
## Normal error learning - for non-binned quality scores
#ATTENTION: Because we decided to run the analyses without the reverse reads, we commented out the lines concerning reverse.

# set the environment
setwd("/Users/anafonseca/Pig_vinny/fastp/forward_reads/")

## ---- SET VARIABLES ----

## modify these variables:

# path to filtered and cleaned reads
CLEANEDPATH = "/Users/anafonseca/Pig_vinny/fastp/forward_reads/"

# path to output tables
OUTPATH = "/Users/anafonseca/Pig_vinny/fastp/forward_reads/output"

# database to silva training set
# downloaded from: https://benjjneb.github.io/dada2/training.html
DB = "/Users/anafonseca/Pig_vinny/silva_data/silva_nr99_v138.1_train_set.fa.gz"

# paired end patterns
FILTEREDF = "_r1_fastp.fq"
#FILTEREDR = "_r2_fastp.fq"

# ---- install and read data ----

### packages must be previously installed
require(dada2)
require(tidyr)
require(phyloseq)
require(magrittr, warn.conflicts = FALSE)

## test that pathway works
#if(!list.files(CLEANEDPATH)) {
#  cat("Can't read file pathway or files are not present")
#}

## ---- core dada algorithm ----

# get forward and reverse reads
forward <- sort(list.files(CLEANEDPATH, pattern = FILTEREDF, full.names = TRUE))
#reverse <- sort(list.files(CLEANEDPATH, pattern = FILTEREDR, full.names = TRUE))

# check to make sure that the lengths of both files are the same and that they match
fwdNames <- sapply(strsplit(basename(forward), FILTEREDF), `[`, 1)
#revNames <- sapply(strsplit(basename(reverse), FILTEREDR), `[`, 1)

# error catch
#if(length(fwdNames) != length(revNames)) {
#  stop("The number of forward and reverse files do not match.")
#} else {

# if(any(!fwdNames%in% revNames)) {

#   stop("Forward and reverse reads are out of order.")
# }
#}

# learn errors for forward and reverse
errF <- learnErrors(
  forward,
  multithread = TRUE,
  verbose = TRUE
)
#errR <- learnErrors(
#  reverse,
# multithread = TRUE,
# verbose = TRUE
#)

# save progress
save.image(file = paste0(OUTPATH, "/error-learning.RData"))

# plot forward errors
pdf(paste0(OUTPATH, "/forward-errorplot.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

# plot reverse errors
#pdf(paste0(OUTPATH, "/reverse-errorplot.pdf"))
#plotErrors(errR, nominalQ = TRUE)
#dev.off()

#derep post-errors

### packages must be previously installed
require(dada2)
require(tidyr)
require(phyloseq)
require(magrittr)


#load post-error/filtering from ben-monotonicity.R
load(paste0(OUTPATH, "/error-learning.RData"))

#dereplicate
derepFs <- derepFastq(forward, verbose=TRUE)
#derepRs <- derepFastq(reverse, verbose=TRUE)

# save image
save.image(file = paste0(OUTPATH, "/derep-image.RData"))

#dada
dadaForward <- dada(derepFs,
                    err=errF,
                    multithread=TRUE,
                    verbose = TRUE)
#dadaReverse <- dada(derepRs,
# err=errR,
# multithread=TRUE,
# verbose = TRUE)

# print finished message
cat("done!")

# save image
save.image(file = paste0(OUTPATH, "/derep-dada2.RData"))

## finishing the dada2 pipeline
# post -dada

## --- CHANGE THESE VARIABLES ----

# set working directory
#OUTPATH <- "/Users/anafonseca/Pig_vinny/fastp/test/output"

# load previous RData
load(paste0(OUTPATH, "/derep-dada2.RData"))

## path to reference database
# downloaded from: https://benjjneb.github.io/dada2/training.html
DB = "/Users/anafonseca/Pig_vinny/silva_data/silva_nr99_v138.1_train_set.fa.gz"

### --- CODE ----

require(dada2)

# merge paired reads
#mergers <- mergePairs(dadaF = dadaForward,
#derepF = forward,
#  dadaR = dadaReverse,
# derepR = reverse,
# verbose = TRUE)



# save intermediate progress
#save.image(paste0(OUTPATH, "/working-image.RData"))

# construct sequence table of ASVs
seqtab <- makeSequenceTable(samples = dadaForward)

# save intermediate progress
save.image(paste0(OUTPATH, "/working-image.RData"))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(unqs = seqtab,
                                    method = "consensus",
                                    verbose = TRUE)

# save intermediate progress
save.image(paste0(OUTPATH, "/working-image.RData"))

# assign taxonomy using the Silva database
tax <- assignTaxonomy(seqs = seqtab.nochim,
                      refFasta = DB,
                      verbose = TRUE)

# save intermediate progress
save.image(paste0(OUTPATH, "/working-image.RData"))

## This is the end of the dada2 algorithm

## ---- WRITE OUTPUT ----

# write ASV table to file
write.table(seqtab.nochim, file = paste0(OUTPATH, "/asv-table.txt"), sep = "\t", row.names = FALSE)

# save ASB table as RData object
save(seqtab.nochim, file = paste0(OUTPATH, "/asv-table.RData"))

# write taxonomy table to file
write.table(tax, file = paste0(OUTPATH, "/tax-table.txt"), sep = "\t", row.names = FALSE)

# save tax table as RData object
save(tax, file = paste0(OUTPATH, "/tax-table.RData"))

# save all objects
save.image(file = paste0(OUTPATH, "/all-dada-objects.RData"))

## ---- track reads through the pipeline ----

# define function (from dada2 tutorial)
getN <- function(x) sum(getUniques(x))

# create dataframe
track <- cbind(sapply(dadaForward, getN), rowSums(seqtab.nochim))

# change names
colnames(track) <- c("denoisedF", "nonchim")
rownames(track) <- fwdNames

# write to file
write.table(track, file = paste0(OUTPATH, "track-reads.txt"), sep = "\t")

# print finished message
cat("done!")
