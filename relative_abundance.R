#load packages
library(phyloseq)
library(microViz)
library(dplyr)
library(stringr)
library(tidylog)
library(tidyverse)
library(data.table)


# set directory
setwd("/Users/anafonseca/OneDrive - The Pennsylvania State University/Pig_vinny/fastp/forward_reads/output/")

# load data
ps <- readRDS("ps-decontam-filtered-counts.rds")

#get otu and taxa tables
counts <- as.data.frame(ps@otu_table)
tax <- as.data.frame(ps@tax_table)

## use the summarize_taxa function
summarize_taxa <- function(features, taxonomy) {
  
  taxlevels<-c("Kingdom","Phylum","Class","Order","Family","Genus")
  
  if(missing(features)){stop("Feature table not provided")}
  if(missing(taxonomy)){stop("taxonomy table not provided")}
  if(sum(colnames(taxonomy) %in% taxlevels)!=6){stop("Taxonomy does not contain expected columns containing Kingdom,Phylum,Class,Order,Family,Genus.")}
  
  output<-list()
  
  for(lvl in taxlevels){
    suppressMessages(
      output[[lvl]]<-
        features %>%
        as.data.frame() %>%
        rownames_to_column("FeatureID") %>%
        gather(-FeatureID, key="SampleID", value="Counts") %>%
        left_join(
          taxonomy %>%
            rownames_to_column("FeatureID") %>%
            unite("Taxon", taxlevels[1:grep(lvl, taxlevels)], sep="; ") %>%
            select(FeatureID, Taxon)
        ) %>%
        group_by(SampleID, Taxon) %>%
        summarize(Counts=sum(Counts)) %>%
        ungroup() %>%
        spread(key=SampleID, value=Counts) %>%
        as.data.frame() %>%
        column_to_rownames("Taxon")
    )
  }
  return(output)
}

taxsummary <- summarize_taxa(counts,tax)

Phylum <- taxsummary$Phylum
Family <- taxsummary$Family
Genus <- taxsummary$Genus

#Split the genus column
Genus$Genus <- str_split_fixed(rownames(Genus),";",6)[,6]
colnames(Genus)
Genus$Genus

Phylum$Phylum <- str_split_fixed(rownames(Phylum),";",6)[,2]
colnames(Phylum)
Phylum$Phylum


Family$Family <- str_split_fixed(rownames(Family),";",6)[,5]
colnames(Family)
Family$Family

#Replace NA to Unclassified
Genus$Genus <- str_replace(Genus$Genus, "NA", "Unclassified")

write_csv(Genus, "RelAbundance_table.csv")


Phylum <- aggregate(Phylum[,sapply(Phylum,is.numeric)],Phylum["Phylum"],sum)
writexl::write_xlsx(Phylum, "Phylum_relabundance.xlsx")

Family <- aggregate(Family[,sapply(Family,is.numeric)],Family["Family"],sum)
writexl::write_xlsx(Family, "Family_relabundance.xlsx")

Genus <- aggregate(Genus[,sapply(Genus,is.numeric)],Genus["Genus"],sum)
writexl::write_xlsx(Genus, "Genus_relabundance.xlsx")

##Three tables were generated as .xlsx files then average, percentage and standard deviation were calculated with Excel functions.