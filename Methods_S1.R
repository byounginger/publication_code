# Supplemental R code for:
# Younginger, Stewart and Ballhorn 2019
# Habitat filtering leads to distinct 
# endophyte communities in ferns over a spatial gradient

library(biom)
library(dplyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(dunn.test)
library(reshape2)
library(data.table)

spatial = 'Table_S2.biom'
mapping = 'Table_S3.txt'
data = import_biom(spatial)
envir = import_qiime_sample_data(mapping)
pseq = merge_phyloseq(data,envir)
rank_names(pseq)
colnames(tax_table(pseq)) =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

plant = merge_samples(pseq, 'Plant')

#write.table(sample_data(plant), file = "spatial_map_combined.txt", sep = "\t")
pl_tb = 'Table_S4.txt' 
pl_tb = import_qiime_sample_data(pl_tb)
otus = otu_table(plant)
tax = tax_table(plant)

# Remerge
plant2 = merge_phyloseq(otus,tax,pl_tb)

plant2_clean <- subset_taxa(plant2, Class != "NA")
plant2_clean <- subset_taxa(plant2_clean, Class != "GS18")
plant2_clean <- subset_taxa(plant2_clean, Genus != "NA")
plant2_clean <- subset_taxa(plant2_clean, Genus != "Candida")

sample_sums(plant2_clean)

set.seed(777)
plant2_rare <- rarefy_even_depth(plant2_clean, sample.size = 1229)

lowSamples = c("ESP_N4")
allSamples <- sample_names(plant2_rare)
allSamples <- allSamples[!(allSamples %in% lowSamples)]
plant3_rare = prune_samples(allSamples, plant2_rare) 

sample_data(plant3_rare)$Site <- factor(sample_data(plant3_rare)$Site, levels = c("Coast ferns","Coast neighbors","Coast Range ferns","Coast Range neighbors","MSH control ferns","MSH control neighbors","MSH impacted ferns","MSH impacted neighbors"))

ESP <- subset_samples(plant3_rare, Sec_Site == "ESP")
SSP <- subset_samples(plant3_rare, Sec_Site =="SSP")
GM <- subset_samples(plant3_rare, Sec_Site == "GM")
MSH <- subset_samples(plant3_rare, Sec_Site == "MSH")

ESP <- prune_taxa(taxa_sums(ESP) > 0, ESP)
SSP <- prune_taxa(taxa_sums(SSP) > 0, SSP)
GM <- prune_taxa(taxa_sums(GM) > 0, GM)
MSH <- prune_taxa(taxa_sums(MSH) > 0, MSH)

### NODF nested analysis ###

# Nested analysis major hypotheses:
# 1. Fungal communities in ferns are nested within or derived from neighboring plants
# 2. More species poor sites are nested within or derived from more species rich sites

full <- as.data.frame(otu_table(plant3_rare))
ESP <- as.data.frame(otu_table(ESP))
SSP <- as.data.frame(otu_table(SSP))
GM <- as.data.frame(otu_table(GM))
MSH <- as.data.frame(otu_table(MSH))

# write.table(full, file = "full.csv", sep = ",")
# write.table(ESP, file = "ESP.csv", sep = ",")
# write.table(SSP, file = "SSP.csv", sep = ",")
# write.table(GM, file = "GM.csv", sep = ",")
# write.table(MSH, file = "MSH.csv", sep = ",")

full2 <- read.csv("Table_S5.csv", header = TRUE, row.names = "SampleID")
ESP2 <- read.csv("Table_S6.csv", header = TRUE, row.names = "SampleID")
SSP2 <- read.csv("Table_S7.csv", header = TRUE, row.names = "SampleID")
GM2 <- read.csv("Table_S8.csv", header = TRUE, row.names = "SampleID")
MSH2 <- read.csv("Table_S9.csv", header = TRUE, row.names = "SampleID")

# Hypothesis 1. Fungal communities in ferns are nested within or derived from 
    # neighboring plants
# Parse out each combination of a target fern and its neighboring plants

## MSH impacted ##

MSHN1 <- MSH2[1,]
MSHN1[2,] <- MSH2[11,]

MSHN2 <- MSH2[2,]
MSHN2[2,] <- MSH2[12,]

MSHN3 <- MSH2[3,]
MSHN3[2,] <- MSH2[13,]

MSHN4 <- MSH2[4,]
MSHN4[2,] <- MSH2[14,]

MSHN5 <- MSH2[5,]
MSHN5[2,] <- MSH2[15,]

MSHN6 <- MSH2[6,]
MSHN6[2,] <- MSH2[16,]

MSHN7 <- MSH2[7,]
MSHN7[2,] <- MSH2[17,]

MSHN8 <- MSH2[8,]
MSHN8[2,] <- MSH2[18,]

MSHN9 <- MSH2[9,]
MSHN9[2,] <- MSH2[19,]

MSHN10 <- MSH2[10,]
MSHN10[2,] <- MSH2[20,]

MSHN11 <- MSH2[10,]
MSHN10[2,] <- MSH2[20,]

# Run the NODF analysis with a null model

oecosimu(MSHN1, nestednodf, "swap", nsimul = 999)
oecosimu(MSHN2, nestednodf, "swap", nsimul = 999)
oecosimu(MSHN3, nestednodf, "swap", nsimul = 999)
oecosimu(MSHN4, nestednodf, "swap", nsimul = 999)
oecosimu(MSHN5, nestednodf, "swap", nsimul = 999)
oecosimu(MSHN6, nestednodf, "swap", nsimul = 999)
oecosimu(MSHN7, nestednodf, "swap", nsimul = 999)
oecosimu(MSHN8, nestednodf, "swap", nsimul = 999)
oecosimu(MSHN9, nestednodf, "swap", nsimul = 999)
oecosimu(MSHN10, nestednodf, "swap", nsimul = 999)

## Coast ##

ESPN1 <- ESP2[1,]
ESPN1[2,] <- ESP2[10,]

ESPN2 <- ESP2[2,]
ESPN2[2,] <- ESP2[11,]

ESPN3 <- ESP2[3,]
ESPN3[2,] <- ESP2[12,]

ESPN4 <- ESP2[4,]
ESPN4[2,] <- ESP2[13,]

ESPN5 <- ESP2[5,]
ESPN5[2,] <- ESP2[15,]

ESPN6 <- ESP2[6,]
ESPN6[2,] <- ESP2[16,]

ESPN7 <- ESP2[7,]
ESPN7[2,] <- ESP2[17,]

ESPN8 <- ESP2[8,]
ESPN8[2,] <- ESP2[18,]

ESPN9 <- ESP2[9,]
ESPN9[2,] <- ESP2[19,]

oecosimu(ESPN1, nestednodf, "swap", nsimul = 999)
oecosimu(ESPN2, nestednodf, "swap", nsimul = 999)
oecosimu(ESPN3, nestednodf, "swap", nsimul = 999)
oecosimu(ESPN4, nestednodf, "swap", nsimul = 999)
oecosimu(ESPN5, nestednodf, "swap", nsimul = 999)
oecosimu(ESPN6, nestednodf, "swap", nsimul = 999)
oecosimu(ESPN7, nestednodf, "swap", nsimul = 999)
oecosimu(ESPN8, nestednodf, "swap", nsimul = 999)
oecosimu(ESPN9, nestednodf, "swap", nsimul = 999)

## Coast Range ##

SSPN1 <- SSP2[1,]
SSPN1[2,] <- SSP2[11,]

SSPN2 <- SSP2[2,]
SSPN2[2,] <- SSP2[12,]

SSPN3 <- SSP2[3,]
SSPN3[2,] <- SSP2[13,]

SSPN4 <- SSP2[4,]
SSPN4[2,] <- SSP2[14,]

SSPN5 <- SSP2[5,]
SSPN5[2,] <- SSP2[15,]

SSPN6 <- SSP2[6,]
SSPN6[2,] <- SSP2[16,]

SSPN7 <- SSP2[7,]
SSPN7[2,] <- SSP2[17,]

SSPN8 <- SSP2[8,]
SSPN8[2,] <- SSP2[18,]

SSPN9 <- SSP2[9,]
SSPN9[2,] <- SSP2[19,]

SSPN10 <- SSP2[10,]
SSPN10[2,] <- SSP2[20,]

oecosimu(SSPN1, nestednodf, "swap", nsimul = 999)
oecosimu(SSPN2, nestednodf, "swap", nsimul = 999)
oecosimu(SSPN3, nestednodf, "swap", nsimul = 999)
oecosimu(SSPN4, nestednodf, "swap", nsimul = 999)
oecosimu(SSPN5, nestednodf, "swap", nsimul = 999)
oecosimu(SSPN6, nestednodf, "swap", nsimul = 999)
oecosimu(SSPN7, nestednodf, "swap", nsimul = 999)
oecosimu(SSPN8, nestednodf, "swap", nsimul = 999)
oecosimu(SSPN9, nestednodf, "swap", nsimul = 999)
oecosimu(SSPN10, nestednodf, "swap", nsimul = 999)

## MSH control ##

GMN1 <- GM2[1,]
GMN1[2,] <- GM2[11,]

GMN2 <- GM2[2,]
GMN2[2,] <- GM2[12,]

GMN3 <- GM2[3,]
GMN3[2,] <- GM2[13,]

GMN4 <- GM2[4,]
GMN4[2,] <- GM2[14,]

GMN5 <- GM2[5,]
GMN5[2,] <- GM2[15,]

GMN6 <- GM2[6,]
GMN6[2,] <- GM2[16,]

GMN7 <- GM2[7,]
GMN7[2,] <- GM2[17,]

GMN8 <- GM2[8,]
GMN8[2,] <- GM2[18,]

GMN9 <- GM2[9,]
GMN9[2,] <- GM2[19,]

GMN10 <- GM2[10,]
GMN10[2,] <- GM2[20,]

oecosimu(GMN1, nestednodf, "swap", nsimul = 999)
oecosimu(GMN2, nestednodf, "swap", nsimul = 999)
oecosimu(GMN3, nestednodf, "swap", nsimul = 999)
oecosimu(GMN4, nestednodf, "swap", nsimul = 999)
oecosimu(GMN5, nestednodf, "swap", nsimul = 999)
oecosimu(GMN6, nestednodf, "swap", nsimul = 999)
oecosimu(GMN7, nestednodf, "swap", nsimul = 999)
oecosimu(GMN8, nestednodf, "swap", nsimul = 999)
oecosimu(GMN9, nestednodf, "swap", nsimul = 999)
oecosimu(GMN10, nestednodf, "swap", nsimul = 999)

# Hypothesis 2. More species poor sites are nested within or derived 
    # from more species rich sites
# Separate ferns and neighbors, pool by site and then run nested analysis

neigh_only <- subset_samples(plant3_rare, Type == "neighbors")
ferns_only <- subset_samples(plant3_rare, Type == "ferns")

neigh_only <- prune_taxa(taxa_sums(neigh_only) > 0, neigh_only)
ferns_only <- prune_taxa(taxa_sums(ferns_only) > 0, ferns_only)

neigh_site <- merge_samples(neigh_only, "Sec_Site")
fern_site <- merge_samples(ferns_only, "Sec_Site")

neigh <- as.data.frame(otu_table(neigh_site))
fern <- as.data.frame(otu_table(fern_site))

oecosimu(neigh, nestednodf, "swap", nsimul = 999)
oecosimu(fern, nestednodf, "swap", nsimul = 999)

### Alpha diversity ###

full <- otu_table(plant3_rare)
Hillfull = renyi(full, scales = c(0,1,2), hill = T)

setDT(Hillfull, keep.rownames = "SampleID")
Hillfull2 <- melt(data = Hillfull, id.vars = "SampleID", measure.vars = c("0","1","2"))
# write.table(Hillfull2, file = "hill_V2.csv", sep = ",")

# # Import the rearranged data:
Hillfull2 <- read.csv(file = "Table_S10.csv", header = TRUE, sep = ",")

# Adjust the factors
Hillfull2$Site <- factor(Hillfull2$Site, levels = c("Coast","Coast neighbors","Coast range",
  "Coast range neighbors","MSH control","MSH control neighbors","MSH impacted","MSH impacted neighbors"))

### Stats on Hill numbers ###

# Separate by Hill number
Hillzero <- dplyr::filter(Hillfull2, q.value == "Hill_0")
Hillone <- dplyr::filter(Hillfull2, q.value == "Hill_1")
Hilltwo <- dplyr::filter(Hillfull2, q.value == "Hill_2")

kruskal.test(Hill.Number ~ Site, data = Hillzero)
dunn.test(x = Hillzero$Hill.Number, g = Hillzero$Site, method = "bh")

kruskal.test(Hill.Number ~ Site, data = Hillone)
dunn.test(x = Hillone$Hill.Number, g = Hillone$Site, method = "bh")

kruskal.test(Hill.Number ~ Site, data = Hilltwo)

dunn.test(x = Hilltwo$Hill.Number, g = Hilltwo$Site, method = "bh")

# Plots

h0 <- ggplot(Hillzero, aes(x=Site, y=Hill.Number, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete("Site", labels = c("Coast","Coast\nneighbors","Coast Range",
  "Coast Range\nneighbors","MSH control","MSH control\nneighbors",
  "MSH\nimpacted","MSH impacted\nneighbors")) + theme_classic() + ylab("Richness\n(OTUs)")

h0 + theme(legend.position = "none", axis.text = element_text(colour = "black"),
  axis.title.x = element_blank(), axis.text.x = element_blank())

h1 <- ggplot(Hillone, aes(x=Site, y=Hill.Number, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete("Site", labels = c("Coast","Coast\nneighbors","Coast Range",
  "Coast Range\nneighbors","MSH control","MSH control\nneighbors",
  "MSH\nimpacted","MSH impacted\nneighbors")) + theme_classic() + ylab("Exp. Shannon\nEntropy")

h1 + theme(legend.position = "none", axis.text = element_text(colour = "black"),
           axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  ylim(0,20)

h2 <- ggplot(Hilltwo, aes(x=Site, y=Hill.Number, fill=Site)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete("Site", labels = c("Coast\nferns","Coast\nneighbors","Coast Range\nferns",
  "Coast Range\nneighbors","MSH control\nferns","MSH control\nneighbors",
  "MSH\nimpacted\ferns","MSH impacted\nneighbors")) + theme_classic() + ylab("Inv. Simpson\nIndex")

h2 + theme(legend.position = "none", axis.text = element_text(colour = "black")) +
  ylim(0, 13)

### Beta diversity ###
# Hypothesis 3. Fungal communities in ferns experience less spatial turnover
#   (i.e. lower beta diversity) than in non-fern hosts

full <- as.data.frame(otu_table(plant3_rare))
sampledf_rare <- data.frame(sample_data(plant3_rare))
full_dist <- vegdist(full, method = "bray", binary = TRUE)
beta2 <- betadisper(full_dist, sampledf_rare$Type)
permutest(beta2)

# Significantly less turnover is observed in the fern hosts

# Make the plot:

beta3 <- as.data.frame(beta2$distances)
beta3$group <- beta2$group
names(beta3)[1] <- "distances"
p <- ggplot(beta3, aes(x = group, y = distances, fill = group)) + geom_boxplot()

p + ylab("Distance to centroid") + xlab("Host type") + 
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_text(color = "black", size = 10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") + scale_x_discrete(labels = c("Ferns", "Neighbors"))

### Rarefaction curves ###

ESP <- subset_samples(plant3_rare, Sec_Site == "ESP")
SSP <- subset_samples(plant3_rare, Sec_Site =="SSP")
GM <- subset_samples(plant3_rare, Sec_Site == "GM")
MSH <- subset_samples(plant3_rare, Sec_Site == "MSH")

ESP <- prune_taxa(taxa_sums(ESP) > 0, ESP)
SSP <- prune_taxa(taxa_sums(SSP) > 0, SSP)
GM <- prune_taxa(taxa_sums(GM) > 0, GM)
MSH <- prune_taxa(taxa_sums(MSH) > 0, MSH)

ESP <- as.data.frame(otu_table(ESP))
SSP <- as.data.frame(otu_table(SSP))
GM <- as.data.frame(otu_table(GM))
MSH <- as.data.frame(otu_table(MSH))

SAC_ESP <- specaccum(ESP, permutations = 999, method = "rarefaction")
SAC_SSP <- specaccum(SSP, permutations = 999, method = "rarefaction")
SAC_GM <- specaccum(GM, permutations = 999, method = "rarefaction")
SAC_MSH <- specaccum(MSH, permutations = 999, method = "rarefaction")

SAC_ESP_df <- data.frame(SAC_ESP$sites, SAC_ESP$richness, SAC_ESP$sd)
SAC_SSP_df <- data.frame(SAC_SSP$sites, SAC_SSP$richness, SAC_SSP$sd)
SAC_GM_df <- data.frame(SAC_GM$sites, SAC_GM$richness, SAC_GM$sd)
SAC_MSH_df <- data.frame(SAC_MSH$sites, SAC_MSH$richness, SAC_MSH$sd)

colids <- c("sites", "richness", "sd")

names(SAC_ESP_df) <- colids
names(SAC_SSP_df) <- colids
names(SAC_GM_df) <- colids
names(SAC_MSH_df) <- colids
SAC_ESP_df$Site <- "ESP"
SAC_SSP_df$Site <- "SSP"
SAC_GM_df$Site <- "GM"
SAC_MSH_df$Site <- "MSH"

combo <- rbind(SAC_ESP_df, SAC_SSP_df, SAC_GM_df, SAC_MSH_df)

combo$Site <- as.factor(combo$Site)

combo$Site <- factor(combo$Site, levels = c("ESP", "SSP", "GM", "MSH"))

p1 <- ggplot(combo, aes(x = sites, y = richness, colour = Site)) + geom_line()

p1 + geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), 
  alpha = 1/5) + theme_classic() + ylab("Richness") + 
  xlab("Plants") + theme(axis.text = element_text(size = 10, color = "black")) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"), 
  labels = c("Coast", "Coast Range", "MSH control", "MSH impacted"))

### Ordination plots ###

# By site
q <- plot_ordination(plant3_rare, ordinate(plant3_rare, "NMDS", distance = "bray", trymax = 100, k = 3), "sites", color = "Site", title = "Stress = 0.16") +
  stat_ellipse() + theme_classic()

q + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12), 
          axis.text = element_text(size = 10, color = "black"), legend.position = "none")

# Ferns vs. neighbors
p <- plot_ordination(plant3_rare, ordinate(plant3_rare, "NMDS", distance = "bray", trymax = 100, k = 3), "sites", color = "Type", title = "Stress = 0.16") +
  stat_ellipse() + theme_classic()

p + theme(legend.text = element_text(size = 12), legend.position = "none", 
          axis.text = element_text(size = 10, color = "black"))

# PERMANOVA w/ Bray-Curtis
bray_rare <- phyloseq::distance(plant3_rare, method = "bray")
sampledf_rare <- data.frame(sample_data(plant3_rare))

# By site
adonis2(bray_rare ~ Site, data = sampledf_rare)
beta <- betadisper(bray_rare, sampledf_rare$Site)
permutest(beta)

# Ferns vs. neighbors 
adonis2(bray_rare ~ Type, data = sampledf_rare)
beta <- betadisper(bray_rare, sampledf_rare$Type)
permutest(beta)

### Taxonomy plot ###

rela = transform_sample_counts(plant3_rare, function(x) x/11050) 
      # The denominator in the function is derived from the sum of all reads 
        # per plant type and site

no_low = filter_taxa(rela, function(x) mean(x) > 1e-3, TRUE)

p <- plot_bar(no_low, "Family", fill = "Genus")
p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + 
  facet_wrap( ~ Site, ncol = 2, nrow = 4) + ylab("Relative abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, color = "black"), 
  axis.text.y = element_text(colour = "black"), 
  panel.background = element_blank(), axis.line = element_line(colour = "black"))
