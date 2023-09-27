# load packages
require(tidyverse)
require(phyloseq)
library(ALDEx2)
library(microViz)
library(writexl)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# set directory
setwd("/Users/anafonseca/OneDrive - The Pennsylvania State University/Pig_vinny/fastp/forward_reads/output/")

###Oral Fluid vs Feces

# load phyloseq of counts
ps <- readRDS("fastp/forward_reads/output/ps-decontam-filtered-pscount.rds")

otu_table(ps)
ps

##select two samples
pscount1 <- subset_samples(ps, Type == "oral fluid" | Type == "feces")
sample_data(pscount1)

# aggregate at genus level for counts
ps1 <- pscount1 %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(pscount1)$Type, taxa_rank = "Genus", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "none")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(pscount1))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# make a table of significant corrected p-values
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
save(sig_aldex2, file = "table_sig_taxa_oralxfeces.RData")
write.table(sig_aldex2, file = "table_sig_taxa_oralxfeces.csv", sep = ",", col.names = TRUE, row.names = FALSE)
write_xlsx(sig_aldex2, "table_sig_taxa_oralxfeces.xlsx")

##Feces vs Pen Floor

##select two samples
pscount2 <- subset_samples(ps, Type == "pen floor swab" | Type == "feces")

sample_data(pscount2)

# aggregate at genus level for counts
ps1 <- pscount2 %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(pscount2)$Type, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "none")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)
#no differences were seen

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(pscount1))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# make a table of significant corrected p-values
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
## no differences

##Oral Fluid vs Pen Floor
sample_data(ps)

##select two samples
pscount3 <- subset_samples(ps, Type == "pen floor swab" | Type == "oral fluid")
sample_data(pscount3)


# aggregate at genus level for counts
ps1 <- pscount3 %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(pscount3)$Type, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "none")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(pscount3))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# make a table of significant corrected p-values
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
save(sig_aldex2, file = "table_sig_taxa_oralxpenfloor.RData")
write.table(sig_aldex2, file = "table_sig_taxa_oralxpenfloor.csv", sep = ",", col.names = TRUE, row.names = FALSE)
write_xlsx(sig_aldex2, "table_sig_taxa_oralxpenfloor.xlsx")


## graph 1

par(mfrow=c(1,2))
plot(sig_aldex2$effect, sig_aldex2$we.ep, log="y", cex=0.7, col=rgb(0,0,1,0.2),
     pch=19, xlab="effect", ylab="wi.eBH", main="Effect size plot")
points(sig_aldex2$effect, sig_aldex2$we.eBH, cex=0.7, col=rgb(1,0,0,0.2),
       pch=19)
abline(h=0.05, lty=2, col="grey")
legend(15,1, legend=c("P value", "BH-adjusted"), pch=19, col=c("blue", "red"))

plot(sig_aldex2$diff.btw, sig_aldex2$we.ep, log="y", cex=0.7, col=rgb(0,0,1,0.2),
     pch=19, xlab="Difference", ylab="P value", main="Volcano plot")
points(sig_aldex2$diff.btw, sig_aldex2$we.eBH, cex=0.7, col=rgb(1,0,0,0.2),
       pch=19)
abline(h=0.05, lty=2, col="grey")


##graph 1
dat1 <- readxl::read_excel("fastp/forward_reads/output/sample_type_table_sig_taxa_oralxfeces.xlsx")

ggplot(data=dat1, aes(x=Genus, y=(effect), col=Sample_type, label=Genus)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  #scale_color_brewer() +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black")


##graph 1
dat2 <- readxl::read_excel("fastp/forward_reads/output/table_sig_taxa_oralxpenfloor.xlsx")

ggplot(data=dat2, aes(x=Genus, y=(effect), col=Sample_Name, label=Genus)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  #scale_color_brewer() +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black")


###Lollipop graph
##Feces vs Oral Fluid
ggplot(dat1, aes(x = effect, y = Genus)) +
  geom_segment(aes(yend = Genus), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(colour = Sample_type)) +
  scale_colour_brewer(palette = "Set1", limits = c("Oral Fluid", "Feces"), guide = "none") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  facet_grid(Sample_type ~ ., scales = "free_y", space = "free_y")

##Oral Fluid vs Pen Floor

###Lollipop graph
ggplot(dat2, aes(x = effect, y = Genus)) +
  geom_segment(aes(yend = Genus), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(colour = Sample_Name)) +
  scale_colour_brewer(palette = "Set1", limits = c("Oral Fluid", "Pen Floor"), guide = "none") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  facet_grid(Sample_Name ~ ., scales = "free_y", space = "free_y")


##cleavand dot plot graph
##Feces vs Oral Fluid
p1 <- ggplot(dat1, aes(effect, Genus)) +
  scale_color_manual(values = c("#028571", "#CC6633")) +
  geom_line(aes(group = Genus)) +
  geom_point(aes(color = Sample_type), size = 4) +
  theme_bw(base_size = 15)

p2 <- ggpar(p1, font.ytickslab = "italic")

p2

##cleavand dot plot graph
##Oral Fluid vs Pen Floor
p3 <- ggplot(dat2, aes(effect, Genus)) +
  scale_color_manual(values = c("#CC6633", "#543005" )) +
  geom_line(aes(group = Genus)) +
  geom_point(aes(color = Sample_Name), size = 4) +
  theme_bw(base_size = 15)

p4 <- ggpar(p3, font.ytickslab = "italic")

p4

ggarrange(p2, p4, common.legend = FALSE, nrow = 2)