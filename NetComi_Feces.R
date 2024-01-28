#Load package
library(phyloseq)
library(NetCoMi)

#Load phyloseq object
ps <- readRDS("/Users/viniciusbuiatte/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/2023-Buiatte-Collab/Scripts/ps-decontam-filtered-pscount.rds")

#Transpose OTU table. ASVs need to be columns 
transp_otu <- t(otu_table(ps))

# Update the phyloseq object with the transposed OTU table
ps1 <- phyloseq(tax_table(ps), sample_data(ps), transp_otu)

#Name the sample type you will run the analysis
Feces <- "feces"

# Filter phyloseq based on sample type
Feces_ps <- subset_samples(ps1, Type == Feces)

#Filter low abundance taxa (<10 counts)
Feces_ps_filt <- prune_taxa(taxa_sums(Feces_ps) > 10, Feces_ps)

#Let's filter out OTUs that are not present in at least 20% of our samples

threshold <- 0.2 * nsamples(Feces_ps)
Feces_ps_3 <- filter_taxa(Feces_ps, function(x) sum(x > 0) >= threshold, TRUE)

#agglomerate at genus level

Feces_genus <- tax_glom(Feces_ps_3, taxrank = "Genus")

taxtab <- as(tax_table(Feces_genus), "matrix")

Feces_genus_renamed <- renameTaxa(Feces_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")

# Network construction and analysis
set.seed(101010)
Feces_genus_net <- netConstruct(Feces_genus_renamed,
                          taxRank = "Genus",
                          measure = "sparcc",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)

set.seed(101010)
props_genus_feces <- netAnalyze(Feces_genus_net, clustMethod = "cluster_louvain")

#Summarize network measures - Degree, betweenness, closeness, eigenvector
set.seed(101010)
summary(props_genus_feces)

colors <- c("#440154ff", "#404688ff", "#287c8eff", "#27ad81ff", "#8fd744ff")

#Plot graph
set.seed (10010)
plot(props_genus_feces, 
     nodeColor = "cluster",
     colorVec = colors,
     nodeSize = "eigenvector",
     #shortenLabels = "intelligent",
     #title1 = "Network analysis at the genus level with sparCC", 
     #showTitle = TRUE,
     cexTitle = 1.0,
     cexNodes = 4.0,
     posCol = "cyan4", 
     negCol = "lightsalmon",
     nodeTransp = 20,
     cexLabels = 0.5,
     rmSingles = TRUE,
     labelScale = FALSE,
     sameClustCol = TRUE,
     borderCol = "black")

# Add the legend for estimated correlation
legend(0.7, 2.0, cex = 1.0, title = "Estimated Correlation:", 
       legend = c("+", "-"), lty = 1, lwd = 3, col = c("cyan4", "lightsalmon"), 
       bty = "n", horiz = TRUE)



