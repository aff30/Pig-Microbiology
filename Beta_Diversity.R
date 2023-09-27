##  Beta Diversity Analysis
## AF - Last Edited April 20-2023

# ---- setup ----
# load packages
library(BiocManager)
library(phyloseq)
library(vegan)
library(patchwork)
library(microViz)
library(ggplot2)
library(dplyr)
library(stringr)
library(pairwiseAdonis)
library(scales)

set.seed(200789)

# check if we have NAs
anyNA(tax_table(ps)[,"Phylum"])
anyNA(tax_table(ps)[,"Phylum"])

# tax fix our phyloseq object
ps <- tax_fix(ps)

ps <- ps %>% tax_fix(unknowns = c("Incertae Sedis"))


### ---- BETA DIVERSITY ----

# clr transform phyloseq objects at Genus level
beta1 <- ps %>%
  tax_transform(trans = "clr", rank = "Genus") %>%
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

#ADONIS test
# farm, sample type, and their interaction were included in the model


vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Type*sample_data(beta1)$Farm, permutations = 999)

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999
#vegan::adonis2(formula = psdist ~ phyloseq::sample_data(beta1)$Type * sample_data(beta1)$Farm, permutations = 999)
#Df SumOfSqs      R2       F Pr(>F)
#phyloseq::sample_data(beta1)$Type                           2    20301 0.23543 22.5999  0.001 ***
#sample_data(beta1)$Farm                                     4    15802 0.18325  8.7955  0.001 ***
#phyloseq::sample_data(beta1)$Type:sample_data(beta1)$Farm   8     7910 0.09172  2.2013  0.001 ***
#Residual                                                   94    42219 0.48960
#Total                                                     108    86231 1.00000

#https://github.com/joey711/phyloseq/issues/201

##  PCA plot
p1 <- ps %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "Type", plot_taxa = 1:5, size = 4, shape = "Farm") +
  scale_color_manual(values = c("#028571", "#CC6633", "#543005")) +
  theme_classic() +
  ggtitle("PCA ") +
  #theme(legend.position = "none") +
  labs(caption = "")
p1

ggsave(filename = "PCA-microViz-samplevsfarm_PLOT1.pdf", dpi = 600)
#https://r-graphics.org/recipe-colors-palette-continuous


p2 <- ps %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "Type", plot_taxa = 1:5, size = 4, shape = "Farm") +
  scale_color_manual(values = c("#028571", "#CC6633", "#543005")) +
  stat_ellipse(aes(group = Type, color = Type)) +
  theme_classic() +
  ggtitle("A") +
  # theme(legend.position = "none") +
  labs(caption = "")
p2

ggsave(filename = "PCA-microViz-farmvssample_plot2.pdf", dpi = 600)

# clr transform phyloseq objects at Genus level
beta1 <- ps %>%
  tax_transform(trans = "clr", rank = "Genus") %>%
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")
vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Type, permutations = 999)

##FARM
# clr transform phyloseq objects at Genus level
beta1 <- ps %>%
  tax_transform(trans = "clr", rank = "Genus") %>%
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")
vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Farm, permutations = 999)


# Pairwise comparison using pairwiseAdonis
data <- phyloseq::sample_data(beta1)
data
pairwiseAdonis::pairwise.adonis(psdist, data$Farm)

#Farm
#pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#1 vs 5  1  4928.312 6.595669 0.13857708   0.001       0.01   *
#1 vs 2  1  5619.112 8.676428 0.17465885   0.001       0.01   *
#1 vs 3  1  6245.276 9.631134 0.19022157   0.001       0.01   *
#1 vs 4  1  3760.217 5.957904 0.12687756   0.001       0.01   *
#5 vs 2  1  3037.856 4.109925 0.08913320   0.002       0.02   .
#5 vs 3  1  2760.020 3.730016 0.08156604   0.003       0.03   .
#5 vs 4  1  3622.767 5.010430 0.10658125   0.002       0.02   .
#2 vs 3  1  3152.394 4.904581 0.10456507   0.001       0.01   *
#2 vs 4  1  2551.518 4.076942 0.08848118   0.002       0.02   .
#3 vs 4  1  3920.980 6.257163 0.12966288   0.001       0.01   *

# Pairwise comparison using pairwiseAdonis
data <- phyloseq::sample_data(beta1)
pairwiseAdonis::pairwise.adonis(psdist, data$Type)

#pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#feces vs oral fluid  1 18372.679 29.316105 0.23208525   0.001      0.003   *
#feces vs pen floor swab  1  2478.977  4.284185 0.06990686   0.001      0.003   *
#oral fluid vs pen floor swab  1  4127.572  6.285506 0.09777486   0.001      0.003   *

### PCA for  farms only

p3 <- ps %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "Type", plot_taxa = 1:5, size = 4) +
  scale_color_manual(values = c("#028571", "#CC6633", "#543005")) +
  stat_ellipse(aes(group = Type, color = Type)) +
  theme_classic() +
  ggtitle("B") +
  # theme(legend.position = "none") +
  labs(caption = "")
p3

ggsave(filename = "PCA-microViz-sample_plot3.pdf", dpi = 600)


p4 <- ps %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "Farm", plot_taxa = 1:5, size = 4, shape = "Farm") +
  scale_color_manual(values = c("#028571", "#dfc27d", "#FFCC33", "#CC6633", "#543005")) +
  #stat_ellipse(aes(group = "Farm", color = "Farm")) +
  theme_classic() +
  ggtitle("C") +
  # theme(legend.position = "none") +
  labs(caption = "")
p4

ggsave(filename = "PCA-microViz-farmellipse_plot3.pdf", dpi = 600)

p5 <- ps %>%
  tax_transform("clr", rank = "Genus") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "Farm", shape = "Farm", size = 4) +
  scale_colour_brewer(palette = "Dark2")

p5

ggsave(filename = "PCA-microViz-farm_noellipse_plot3.pdf", dpi = 600)

###_______________________________________________________________