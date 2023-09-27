# Decontamination and filtering
#adapted by ana 07/20/22

# load packages
require(tidyverse)
require(phyloseq)
require(dada2)
# install latest version of microViz
#devtools::install_github("david-barnett/microViz@0.7.1")
require(microViz)
#require(BiocManager)
#BiocManager::install("decontam")
require(decontam)
require(ggpubr)

setwd("/Users/anafonseca/OneDrive - The Pennsylvania State University/Pig_vinny/fastp/forward_reads/output/")
#
# get data
load("ps-raw.rds")

## ---- basic info ----
ps <- psraw

sample_data(ps)

ntaxa(ps) # 8058 taxa

get_taxa_unique(ps, "Phylum") # 24 different phyla including NA
get_taxa_unique(ps, "Kingdom") # 3 kingdoms

# how many controls at each step
ps %>% samdat_tbl() %>% group_by(Type) %>% summarize(n())


# ---- filter for Bacteria & remove NA phyla ----

#look for mito or chloroplasts
sort(get_taxa_unique(ps, "Order"))
sort(get_taxa_unique(ps, "Order"))

#remove them if there
psfilt <- subset_taxa(ps, Order != "Chloroplast")
psfilt <- subset_taxa(psfilt, Family != "Mitochondria")

# get only bacteria
psf <- subset_taxa(psfilt, Kingdom == "Bacteria")

# remove NA phylum
psf1 <- subset_taxa(psf, !is.na(Phylum) & !Phylum %in% c("", "NA"))

# validate
psf2 <- tax_fix(psf1)
psf3 <- phyloseq_validate(psf2, remove_undetected = TRUE)

get_taxa_unique(psf3, "Kingdom")
otu_table(psf3)

pscount <- psf3
pscount
otu_table(pscount)

saveRDS(pscount, file = "ps-decontam-filtered-pscount.rds")

## ---- filter for relative abundance ----
# transform to relative abundance
pst <- transform_sample_counts(pscount, function(x) x / sum(x) )


# remove taxa with total relative abundance less than 10e-5
psr <- filter_taxa(pst, function(x) mean(x) > 1e-5, TRUE)


## ---- get phyloseq of counts ----
## we need the count data with the taxa present in the decontaminated and filtered PS
psrelab <- subset_taxa(psr, taxa_names(psr) %in% taxa_names(psr))

psrelab

# save
saveRDS(psrelab, file = "ps-decontam-filtered-relabund.rds")
otu_table(psrelab)

## ---- Relative abundance barplot ----
## all samples
psrelab %>%
  tax_fix() %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 10,
    merge_other = TRUE,
  )+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_wrap(~ Type)

otu_table(psrelab)
save.image(file = "working_decontam.RData")