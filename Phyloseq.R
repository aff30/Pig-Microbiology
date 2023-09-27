## Creating a phyloseq object

# Ana 07/20/22

# install (if necessary) and load packages
require(tidyverse)
#install phyloseq
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")
require(phyloseq)
#install DADA2
#install.packages("devtools")
library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16")
require(dada2)
# install latest version of microViz
#devtools::install_github("david-barnett/microViz@0.7.1")
require(microViz)
require(BiocManager)
#BiocManager::install("decontam")
require(decontam)
require(ggpubr)

# set directory
setwd("/Users/anafonseca/OneDrive - The Pennsylvania State University/Pig_vinny/fastp/forward_reads/output/")

### ---- phyloseq ----
# LOAD ASV table
load("asv-table.RData")

# load tax table
load("tax-table.RData")

# matrices in R
dim(seqtab.nochim)
dim(tax)

# create phyloseq object
ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows = TRUE),
               tax_table(tax))

# look at the phyloseq object
ps

# get number of taxa
ntaxa(ps)

#get taxa ranks
rank_names(ps)

# access the data "slots" with @
head(ps@tax_table)
head(ps@otu_table)

# fix ASV names
### from dada2 tutorial: fix ASV names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
DNAps <- merge_phyloseq(ps, dna)
taxa_names(DNAps) <- paste0("ASV", seq(ntaxa(DNAps)))
DNAps

# get taxa names
head(taxa_names(DNAps))
head(sample_names(DNAps))


## ---- metadata ----
# read in metadata file

dat <- read.delim("/Users/anafonseca/OneDrive - The Pennsylvania State University/Pig_vinny/metadata/Metadata.txt", header = TRUE)

# look at data
head(dat)
str(dat)

# fix sample names to get ONLY the sample ID
#names <- sapply(str_split(sample_names(ps), "_"), `[`, 1) # I had an error with telling me about duplicates
names <- str_remove_all(sample_names(DNAps), "_S(\\d){1,3}_r1_fastp.fq")


# change sample names to NAMES
sample_names(DNAps) <- names

# format our data to add to phyloseq
sampdf <- dat %>%
  column_to_rownames(var = "Project_ID")


# add to phyloseq
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows= FALSE),
               sample_data(sampdf),
               tax_table(tax))

## build phyloseq
sample_data(DNAps) <- sampdf
DNAps


psraw <- DNAps

# save as Rimage
save(psraw, file = "ps-raw.rds")

psraw
otu_table(psraw)