#Alpha Diversity Analysis
# April 20, 2023 - AF

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# set directory
setwd("/Users/anafonseca/OneDrive - The Pennsylvania State University/Pig_vinny/fastp/forward_reads/output/")

# load data
ps <- readRDS("ps-decontam-filtered-counts.rds")

# check phyloseq
ps


# create data frames ----
adivA <- data.frame(
  #"Observed" = phyloseq::estimate_richness(ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(ps, measures = "Simpson"),
  "Type" = phyloseq::sample_data(ps)$Type,
  "Farm" = phyloseq::sample_data(ps)$Farm
)

testinteraction <- aov(Shannon ~ Type*Farm, data = adivA)
summary(testinteraction)

# repeated measures anova ----
# https://m-clark.github.io/docs/mixedModels/anovamixed.html

anovaA <- aov(Shannon ~ Type*Farm, data = adivA)
summary(anovaA)

#Df Sum Sq Mean Sq F value   Pr(>F)
#Type          2 15.115   7.558  40.150 1.28e-13 ***
# Farm          1  0.175   0.175   0.931    0.337
#Type:Farm     2  0.385   0.192   1.022    0.364
#Residuals   103 19.388   0.188

# Farm did not come up significant so I checked
testfarm <- aov(Shannon ~ Farm, data = adivA)
summary(testfarm) # 0.521

testSample <- aov(Shannon ~ Type, data = adivA)
summary(testSample) # 1.04e-13 ***

TukeyHSD(testSample)

#Fit: aov(formula = Shannon ~ Type, data = adivA)

#$Type
#diff        lwr        upr     p adj
#oral fluid-feces          -0.7714864 -0.9787736 -0.5641992 0.0000000
#pen floor swab-feces      -0.1836085 -0.5414295  0.1742126 0.4441999
#pen floor swab-oral fluid  0.5878780  0.2306639  0.9450920 0.0004722


## box plots -----
sample_data(ps)$"Farm" <- factor(sample_data(ps)$"Farm",
                                 levels = c("1", "2", "3", "4", "5"))
sample_data(ps)$"Type" <- factor(sample_data(ps)$"Type",
                                 levels = c("oral fluid", "feces", "pen floor swab"))

A <- plot_richness(ps, x="Farm", measures=c("Shannon"), color = "Type") +
  geom_boxplot() +
  scale_color_manual(values = c("#028571", "#CC6633", "#543005")) +
  theme_classic() +
  theme(legend.position = "none")


A



ggsave(filename = "alpha-sampletypevsfarm_.pdf", dpi = 600)
ggsave(filename = "alpha-sampletype_.pdf", dpi = 600)

# lineplot ----

adivA %>%
  ggline(x = "Farm", y = "Shannon", group = "Type", color = "Type", palette = c("#028571", "#CC6633", "#543005", "#dfc27d", "#FFCC33"), add = "mean_sd", error.plot = "errorbar")

ggsave(filename = "geomline-plot-sampletype*farm.pdf", dpi = 600)

data_new <- adivA # Duplicate data
data_new$Type <- factor(data_new$Type, # Reorder factor levels
                        c("feces", "pen floor swab", "oral fluid"))


#boxplot
B <- ggboxplot(data_new, x = "Type", y = "Shannon",
               add = "jitter", fill = "Type", palette = c("#028571","#543005", "#CC6633"), ylim=c(3,6))
print(B)

ggsave(filename = "alpha_boxplot-sampletype.pdf", dpi = 600)

## select one type and then do the graphs
#oral fluid

ps1 <- ps %>%
  ps_filter(
    Type == "oral fluid"
  )
ps1

# create data frames ----
adivA1 <- data.frame(
  #"Observed" = phyloseq::estimate_richness(ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps1, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(ps1, measures = "Simpson"),
  "Type" = phyloseq::sample_data(ps1)$Type,
  "Farm" = phyloseq::sample_data(ps1)$Farm
)


testfarm_OF <- aov(Shannon ~ Farm, data = adivA1)
summary(testfarm_OF) # 0.521

TukeyHSD(testfarm_OF)

## box plots -----
sample_data(ps1)$"Farm" <- factor(sample_data(ps1)$"Farm",
                                  levels = c("1", "2", "3", "4", "5"))
#sample_data(ps1)$"Type" <- factor(sample_data(ps)$"Type",
#levels = c("oral fluid", "feces", "pen floor swab"))

A1 <- plot_richness(ps1, x="Farm", measures=c("Shannon"), title = "Oral Fluid", color  = "Farm") +
  geom_boxplot() +
  scale_color_manual(values = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600")) +
  theme_classic() +
  theme(legend.position = "none")


A1

#boxplot
box1 <- ggboxplot(adivA1, x = "Farm", y = "Shannon",
                  add = "jitter", fill = "Farm", palette = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600"), title = "Oral Fluid", ylim=c(3,6))
print(box1)

ggsave(filename = "boxplot1_oralfluidxfarm.pdf", dpi = 600)

#violin
p1 <- ggviolin(adivA1, x = "Farm", y = "Shannon",
               add = "boxplot", fill = "Farm", palette = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600"), title = "Oral Fluid")
print(p1)

ggsave(filename = "violin1_oralfluidxfarm.pdf", dpi = 600)

# feces

ps2 <- ps %>%
  ps_filter(
    Type == "feces"
  )
ps2

# create data frames ----
adivA2 <- data.frame(
  #"Observed" = phyloseq::estimate_richness(ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps2, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(ps2, measures = "Simpson"),
  "Type" = phyloseq::sample_data(ps2)$Type,
  "Farm" = phyloseq::sample_data(ps2)$Farm
)


testfarm_FE <- aov(Shannon ~ Farm, data = adivA2)
summary(testfarm_FE) # 0.521

TukeyHSD(testfarm_FE)

## box plots -----
sample_data(ps2)$"Farm" <- factor(sample_data(ps2)$"Farm",
                                  levels = c("1", "2", "3", "4", "5"))
#sample_data(ps1)$"Type" <- factor(sample_data(ps)$"Type",
#levels = c("oral fluid", "feces", "pen floor swab"))

A2 <- plot_richness(ps2, x="Farm", measures=c("Shannon"), title = "Feces", color  = "Farm") +
  geom_boxplot() +
  scale_color_manual(values = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600")) +
  theme_classic() +
  theme(legend.position = "none")


A2

#boxplot
box2 <- ggboxplot(adivA2, x = "Farm", y = "Shannon",
                  add = "jitter", fill = "Farm", palette = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600"), title = "Feces", ylim=c(3,6))
print(box2)

ggsave(filename = "boxplot2_fecesxfarm.pdf", dpi = 600)

#violin
p2 <- ggviolin(adivA2, x = "Farm", y = "Shannon",
               add = "boxplot", fill = "Farm", palette = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600"), title = "Feces")
print(p2)

ggsave(filename = "violin2_fecesxfarm.pdf", dpi = 600)

# Floor
ps3 <- ps %>%
  ps_filter(
    Type == "pen floor swab"
  )
ps3

# create data frames ----
adivA3 <- data.frame(
  #"Observed" = phyloseq::estimate_richness(ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps3, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(ps3, measures = "Simpson"),
  "Type" = phyloseq::sample_data(ps3)$Type,
  "Farm" = phyloseq::sample_data(ps3)$Farm
)

# Farm did not come up significant so I checked
testfarm_PF <- aov(Shannon ~ Farm, data = adivA3)
summary(testfarm_PF) # 0.521

TukeyHSD(testfarm_PF)

## box plots -----
sample_data(ps3)$"Farm" <- factor(sample_data(ps3)$"Farm",
                                  levels = c("1", "2", "3", "4", "5"))
#sample_data(ps1)$"Type" <- factor(sample_data(ps)$"Type",
#levels = c("oral fluid", "feces", "pen floor swab"))

A3 <- plot_richness(ps3, x="Farm", measures=c("Shannon"), title = "Floor", color  = "Farm") +
  geom_boxplot() +
  scale_color_manual(values = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600")) +
  theme_classic() +
  theme(legend.position = "none")


A3

#boxplot
box3 <- ggboxplot(adivA3, x = "Farm", y = "Shannon",
                  add = "jitter", fill = "Farm", palette = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600"), title = "Floor", ylim=c(3,6))
print(box3)

ggsave(filename = "boxplot3_floorxfarm.pdf", dpi = 600)

#violin
p3 <- ggviolin(adivA3, x = "Farm", y = "Shannon",
               add = "boxplot", fill = "Farm", palette = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600"), title = "Floor")
print(p3)

ggsave(filename = "violin3_floorxfarm.pdf", dpi = 600)

#anova

anovaA1 <- aov(Shannon ~ Type, data = adivA)
summary(anovaA1)

anova2 <- aov(Shannon ~ Farm, data = adivA)
summary(anova2)