# set directory
setwd("/Users/anafonseca/OneDrive - The Pennsylvania State University/Pig_vinny/fastp/forward_reads/output/")

# load data
load("working_decontam.RData")


library(phylosmith)
options(width = 100)
library(microViz)
library(phyloch)
library(phyloseq)
library(ggplot2)
library(patchwork) # for arranging groups of plots
knitr::opts_chunk$set(fig.height = 6, fig.width = 9)

hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Family"

# Sort phyloseq at lower, and then higher ranks
ps <- ps %>% tax_fix() %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" families will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(ps)[, hueRank]),
  shade = as.vector(tt_get(ps)[, shadeRank]),
  counts = taxa_sums(otu_get(ps))
)

hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )


hierarchicalPalMatrix <- matrix(
  data = sapply(
    X = seq(from = 30, to = 75, length.out = nShades),
    FUN = function(l) scales::hue_pal(l = l, h.start = 30)(n = nHues)
  ),
  byrow = TRUE, ncol = nHues
)
hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))

hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))

tax_palette_plot(hierarchicalPal) +
  theme(axis.text.y.left = element_text(family = "mono"))


ps %>%
  ps_get() %>%
  tax_mutate("Phylum: Family" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  ps_select(Farm, Type) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Type") %>%
  comp_barplot(
    tax_level = "Phylum: Family", n_taxa = length(hierarchicalPal),
    tax_order = "asis", palette = hierarchicalPal, bar_width = 0.975
  ) +
  coord_flip() +
  theme(legend.text = element_text(family = "mono")) # for text alignment


ggsave("new_plot_phylum.pdf")

### compoositional

pscount %>% tax_fix() %>%
  ps_select(Type, Farm) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Type") %>%
  comp_barplot(tax_level = "Phylum", n_taxa = 12, bar_width = 0.8) +
  coord_flip() + labs(x = NULL, y = NULL)

ggsave("plot_phylum2.pdf")