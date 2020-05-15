## Processing of the temporal MS

setwd('/Users/brett/Documents/Temporal_turnover/Informatics')

# library(biom)
# library(dplyr)
# library(ggbiplot)
# library(MASS)
library(phyloseq)
# library(qiimer)
library(vegan)
# library(gridExtra)
library(ggplot2)
library(dplyr)

##################

### Figure height-width variables for Ecology (ESA)
prop = 5/7      
s_w = 3.00
s_h = s_w*prop
h_w = 4.50
h_h = h_w*prop
d_w = 6.00
d_h = d_w*prop
### Figure height-width variables for Ecology (ESA)

#1. Load files

# OTU analysis:
# otu <- 'merged_otutab_V4.txt'
# map <- '/Users/brett/Documents/Temporal_turnover/Informatics/map_file16.txt'
# tax <- 'merged_otus_V3.sintax'

# ZOTU (ASV or ESV) analysis:
otu <- 'merged_zotutab_V3.txt'
map <- 'map_file16.txt'
tax <- 'merged_zotus_V4_rdp.sintax'

full_otu <- read.delim(otu, header = TRUE, row.names = 1)
full_tax <- read.delim(tax, header = FALSE, row.names = 1)
tax_mat <- as.matrix(full_tax)

full_otu2 <- otu_table(full_otu, taxa_are_rows = TRUE)
full_tax2 <- tax_table(tax_mat)
full_map <- import_qiime_sample_data(map)

pseq <- merge_phyloseq(full_otu2, full_tax2, full_map)

apr_may <- subset_samples(pseq, Month == 'April' | Month == 'May')

jun_jul <- subset_samples(pseq, Month == 'June' | Month == 'July')

aug_sep <- subset_samples(pseq, Month == 'August' | Month == 'September')

oct_nov <- subset_samples(pseq, Month == 'October' | Month == 'November')

dec_jan <- subset_samples(pseq, Month == 'December' | Month == 'January')

# Now check sample sizes and make sure the controls are in respective months
sample_data(pseq)
# 795 samples

sample_data(apr_may)
# 160 samples inc. one control

sample_data(jun_jul)
# 162 samples inc. two controls

sample_data(aug_sep)
# 157 samples inc. two controls

sample_data(oct_nov)
# 156 samples inc. three controls

sample_data(dec_jan)
# 160 samples inc. four controls

160 + 162 + 157 + 156 + 160

# All samples are accounted for. Now extract the otu tables and sum controls,
# then substract from rest of samples

apr_may_otu <- otu_table(apr_may)
am_df <- as.data.frame(apr_may_otu)
am_minus <- am_df[,1:159] - am_df[,160]
am_minus[,1:159][am_minus[,1:159] < 0] <- 0
am_clean <- otu_table(am_minus, taxa_are_rows = TRUE)
am_sam <- subset_samples(sample_data(apr_may), Plant != 'control')
am_tax <- tax_table(apr_may)
am_new <- merge_phyloseq(am_sam, am_tax, am_clean)

##
jun_jul_otu <- otu_table(jun_jul)
jj_df <- as.data.frame(jun_jul_otu)
jj_df['all_cont'] <- jj_df['jj.ex.cont'] + jj_df['jj.mock.comm']
jj_minus <- jj_df[,1:160] - jj_df[,163]
jj_minus2 <- jj_minus[,1:160]
jj_minus2[,1:160][jj_minus2[,1:160] < 0] <- 0
jj_clean <- otu_table(jj_minus2, taxa_are_rows = TRUE)
jj_sam <- subset_samples(sample_data(jun_jul), Plant != 'control')
jj_tax <- tax_table(jun_jul)
jj_new <- merge_phyloseq(jj_sam, jj_tax, jj_clean)

##
aug_sep_df <- as.data.frame(otu_table(aug_sep))
aug_sep_df['all_cont'] <- aug_sep_df['as.mock.comm'] + aug_sep_df['as.pcr.cont']
as_minus <- aug_sep_df[,1:155] - aug_sep_df[,158]
as_minus2 <- as_minus[,1:155]
as_minus2[,1:155][as_minus2[,1:155] < 0] <- 0
as_clean <- otu_table(as_minus2, taxa_are_rows = TRUE)
as_sam <- subset_samples(sample_data(aug_sep), Plant != 'control')
as_tax <- tax_table(aug_sep)
as_new <- merge_phyloseq(as_clean, as_sam, as_tax)

##
on_df <- as.data.frame(otu_table(oct_nov))
on_df['all_cont'] <- on_df['on.ex.cont.1'] + on_df['on.ex.cont.2'] + on_df['on.mock.comm']
on_minus <- on_df[,1:153] - on_df[,157]
on_minus2 <- on_minus[,1:153]
on_minus2[,1:153][on_minus2[,1:153] < 0] <- 0
on_clean <- otu_table(on_minus2, taxa_are_rows = TRUE)
on_tax <- tax_table(oct_nov)
on_sam <- subset_samples(sample_data(oct_nov), Plant != 'control')
on_new <- merge_phyloseq(on_clean, on_tax, on_sam)

## 
dj_df <- as.data.frame(otu_table(dec_jan))
dj_df['all_cont'] <- dj_df[,157] + dj_df[,158] + dj_df[,159] + dj_df[,160]
dj_minus <- dj_df[,1:156] - dj_df[,161]
dj_minus[,1:156][dj_minus[,1:156] < 0] <- 0
dj_clean <- otu_table(dj_minus, taxa_are_rows = TRUE)
dj_tax <- tax_table(dec_jan)
dj_sam <- subset_samples(sample_data(dec_jan), Plant != 'control')
dj_new <- merge_phyloseq(dj_clean, dj_tax, dj_sam)

pseq_clean <- merge_phyloseq(am_new, jj_new, as_new, on_new, dj_new)

sort(sample_sums(pseq_clean))

rank_names(pseq_clean)
colnames(tax_table(pseq_clean)) =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Maybe provide a metric on dissimilarity between samples within plant by month
# but merge samples by plant for now
# Can be more strategic about merging these now to maintain sample_data integrity

pseq_merge <- merge_samples(pseq_clean, 'Plant')

sam_mer <- as.data.frame(sample_data(pseq_merge))
sam_mer_rows <- rownames(sam_mer)
sam_mer['New_SampleID'] <- sam_mer_rows

##
sam_mer['New_Month'] <- sub('\\d+\\.(\\w+)', '\\1', sam_mer$New_SampleID)

sam_mer['New_Plant'] <- as.factor(sub('(\\d+)\\.\\w+', '\\1', sam_mer$New_SampleID))

sam_mer <- sam_mer[,5:7]
new_names <- c('SampleID', 'Month', 'Plant')
colnames(sam_mer) <- new_names

mer_tax <- tax_table(pseq_merge)
mer_otu <- otu_table(pseq_merge)

pseq_merge2 <- merge_phyloseq(mer_otu, mer_tax, sam_mer)
head(sample_data(pseq_merge2))

sort(sample_sums(pseq_merge2))

## Remove some spurious OTUs, potential contaminants? 

pseq_merge3 <- prune_taxa(taxa_sums(pseq_merge2) > 1, pseq_merge2) # Maybe adjust this later? 
# pseq_merge3 <- subset_taxa(pseq_merge3, Genus != 'Candida') # Going to hold off on this for now

sort(sample_sums(pseq_merge3))

## Good for now. Will come back and clean more after visualizing some 
# taxonomic plots
set.seed(777)
pseq_rare <- rarefy_even_depth(pseq_merge3, sample.size = 15000, replace = FALSE)
# Lost 11 samples and 134 OTUs

# Plot some rarefaction curves

#APR <- subset_samples(pseq_rare, Month == 'Apr')
MAY <- subset_samples(pseq_rare, Month == 'May')
JUN <- subset_samples(pseq_rare, Month == 'Jun')
JUL <- subset_samples(pseq_rare, Month == 'Jul')
AUG <- subset_samples(pseq_rare, Month == 'Aug')
SEP <- subset_samples(pseq_rare, Month == 'Sep')
OCT <- subset_samples(pseq_rare, Month == 'Oct')
NOV <- subset_samples(pseq_rare, Month == 'Nov')
DEC <- subset_samples(pseq_rare, Month == 'Dec')
JAN <- subset_samples(pseq_rare, Month == 'Jan')

#APR <- prune_taxa(taxa_sums(APR) > 0, APR)
MAY <- prune_taxa(taxa_sums(MAY) > 0, MAY)
JUN <- prune_taxa(taxa_sums(JUN) > 0, JUN)
JUL <- prune_taxa(taxa_sums(JUL) > 0, JUL)
AUG <- prune_taxa(taxa_sums(AUG) > 0, AUG)
SEP <- prune_taxa(taxa_sums(SEP) > 0, SEP)
OCT <- prune_taxa(taxa_sums(OCT) > 0, OCT)
NOV <- prune_taxa(taxa_sums(NOV) > 0, NOV)
DEC <- prune_taxa(taxa_sums(DEC) > 0, DEC)
JAN <- prune_taxa(taxa_sums(JAN) > 0, JAN)

#Apr_otu <- as.data.frame(otu_table(APR))
May_otu <- as.data.frame(otu_table(MAY))
Jun_otu <- as.data.frame(otu_table(JUN))
Jul_otu <- as.data.frame(otu_table(JUL))
Aug_otu <- as.data.frame(otu_table(AUG))
Sep_otu <- as.data.frame(otu_table(SEP))
Oct_otu <- as.data.frame(otu_table(OCT))
Nov_otu <- as.data.frame(otu_table(NOV))
Dec_otu <- as.data.frame(otu_table(DEC))
Jan_otu <- as.data.frame(otu_table(JAN))

#Apr_otu_sac <- specaccum(Apr_otu, permutations = 999, method = 'rarefaction')
May_otu_sac <- specaccum(May_otu, permutations = 999, method = 'rarefaction')
Jun_otu_sac <- specaccum(Jun_otu, permutations = 999, method = 'rarefaction')
Jul_otu_sac <- specaccum(Jul_otu, permutations = 999, method = 'rarefaction')
Aug_otu_sac <- specaccum(Aug_otu, permutations = 999, method = 'rarefaction')
Sep_otu_sac <- specaccum(Sep_otu, permutations = 999, method = 'rarefaction')
Oct_otu_sac <- specaccum(Oct_otu, permutations = 999, method = 'rarefaction')
Nov_otu_sac <- specaccum(Nov_otu, permutations = 999, method = 'rarefaction')
Dec_otu_sac <- specaccum(Dec_otu, permutations = 999, method = 'rarefaction')
Jan_otu_sac <- specaccum(Jan_otu, permutations = 999, method = 'rarefaction')

#Apr_otu_sac_df <- data.frame(Apr_otu_sac$sites, Apr_otu_sac$richness, Apr_otu_sac$sd)
May_otu_sac_df <- data.frame(May_otu_sac$sites, May_otu_sac$richness, May_otu_sac$sd)
Jun_otu_sac_df <- data.frame(Jun_otu_sac$sites, Jun_otu_sac$richness, Jun_otu_sac$sd)
Jul_otu_sac_df <- data.frame(Jul_otu_sac$sites, Jul_otu_sac$richness, Jul_otu_sac$sd)
Aug_otu_sac_df <- data.frame(Aug_otu_sac$sites, Aug_otu_sac$richness, Aug_otu_sac$sd)
Sep_otu_sac_df <- data.frame(Sep_otu_sac$sites, Sep_otu_sac$richness, Sep_otu_sac$sd)
Oct_otu_sac_df <- data.frame(Oct_otu_sac$sites, Oct_otu_sac$richness, Oct_otu_sac$sd)
Nov_otu_sac_df <- data.frame(Nov_otu_sac$sites, Nov_otu_sac$richness, Nov_otu_sac$sd)
Dec_otu_sac_df <- data.frame(Dec_otu_sac$sites, Dec_otu_sac$richness, Dec_otu_sac$sd)
Jan_otu_sac_df <- data.frame(Jan_otu_sac$sites, Jan_otu_sac$richness, Jan_otu_sac$sd)

colids <- c('plant', 'richness', 'sd')

#colnames(Apr_otu_sac_df) <- colids
colnames(May_otu_sac_df) <- colids
colnames(Jun_otu_sac_df) <- colids
colnames(Jul_otu_sac_df) <- colids
colnames(Aug_otu_sac_df) <- colids
colnames(Sep_otu_sac_df) <- colids
colnames(Oct_otu_sac_df) <- colids
colnames(Nov_otu_sac_df) <- colids
colnames(Dec_otu_sac_df) <- colids
colnames(Jan_otu_sac_df) <- colids

#Apr_otu_sac_df$month <- 'Apr'
May_otu_sac_df$month <- 'May'
Jun_otu_sac_df$month <- 'Jun'
Jul_otu_sac_df$month <- 'Jul'
Aug_otu_sac_df$month <- 'Aug'
Sep_otu_sac_df$month <- 'Sep'
Oct_otu_sac_df$month <- 'Oct'
Nov_otu_sac_df$month <- 'Nov'
Dec_otu_sac_df$month <- 'Dec'
Jan_otu_sac_df$month <- 'Jan'

#full_curve <- rbind(Apr_otu_sac_df, May_otu_sac_df, Jun_otu_sac_df, Jul_otu_sac_df, Aug_otu_sac_df, Sep_otu_sac_df, Oct_otu_sac_df, Nov_otu_sac_df, Dec_otu_sac_df, Jan_otu_sac_df)
full_curve <- rbind(May_otu_sac_df, Jun_otu_sac_df, Jul_otu_sac_df, Aug_otu_sac_df, Sep_otu_sac_df, Oct_otu_sac_df, Nov_otu_sac_df, Dec_otu_sac_df, Jan_otu_sac_df)

full_curve$month <- factor(full_curve$month, levels = c('May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'))

p1 <- ggplot(full_curve, aes(x = plant, y = richness, colour = month)) + geom_line() + xlim(1,20)
p1 + theme_classic() + xlab('Plants') + ylab('Richness') +
  theme(legend.key.height = unit(0.1, 'in')) + scale_color_discrete('Months')

ggsave('../R/Figures/rare_v1.pdf', width = s_w, height = s_h, units = 'in')

# Think instead I should set the xlim to 1 # Will work on this further at a later point
# Also, April is really the outlier due to low sample sizes.
  # I think I should be removing it from the analysis, but will think about it some more

# Ordination

pseq_rare_noApr <- subset_samples(pseq_rare, Month != 'Apr')

sample_data(pseq_rare_noApr)$Month <- factor(sample_data(pseq_rare_noApr)$Month, 
                                       levels = c('May', 'Jun', 'Jul','Aug', 'Sep', 
                                                  'Oct', 'Nov', 'Dec', 'Jan'))

# sample_data(pseq_rare)$Month <- factor(sample_data(pseq_rare)$Month, 
#                                        levels = c('Apr', 'May', 'Jun', 'Jul','Aug', 'Sep', 
#                                                   'Oct', 'Nov', 'Dec', 'Jan'))

p2 <- plot_ordination(pseq_rare_noApr, ordinate(pseq_rare_noApr, 'NMDS', distance = 'bray', trymax = 100), 
                      'sites', color = 'Month', title = 'Stress = 0.15') +
  stat_ellipse() + theme_classic()

p2

ggsave('../R/Figures/ord_v1.pdf', width = d_w, height = d_h, units = 'in')

## Taxonomy

rela = transform_sample_counts(pseq_rare, function(x) x/300000)

no_low = filter_taxa(rela, function(x) mean(x) > 5e-5, TRUE)

p3 <- plot_bar(no_low, "Family", fill = "Genus")
p3 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + 
  facet_wrap( ~ Month, ncol = 5, nrow = 2) + ylab("Relative abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Will need to revise the tax table to reflect the new taxonomy of the genus
# and other taxonomic classifications

# What else, alpha diversity plots
# Mock community analysis (do this after I decide on whether to use OTUs or 
# ASVs) # Maybe do a PERMANOVA or DISSIMILARITY test between the two? 

## 2020-03-13

# Processing the MC samples:

controls <- subset_samples(pseq, Plant == 'control')
nsamples(controls)
sample_data(controls)
colnames(otu_table(controls))

# Substract the extraction and control seqs from each run
otu_cont <- as.data.frame(otu_table(controls))
colnames(otu_cont)
otu_cont['as.mock.minus'] <- otu_cont['as.mock.comm'] - otu_cont['as.pcr.cont']
otu_cont['jj.mock.minus'] <- otu_cont['jj.mock.comm'] - otu_cont['jj.ex.cont']
otu_cont['on.ex.sum'] <- otu_cont['on.ex.cont.1'] + otu_cont['on.ex.cont.2']
otu_cont['on.mock.minus'] <- otu_cont['on.mock.comm'] - otu_cont['on.ex.sum']
otu_cont['dj.ex.sum'] <- otu_cont['dj.ex.cont.2'] + otu_cont['dj.ex.cont.3'] + 
  otu_cont['dj.pcr.cont']
otu_cont['dj.mock.minus'] <- otu_cont['dj.mock.comm'] - otu_cont['dj.ex.sum']

otu_cont2 <- otu_cont[,13:18]
colnames(otu_cont2)
keepers <- c("as.mock.minus", "jj.mock.minus", "on.mock.minus", "dj.mock.minus")
otu_cont2 <- otu_cont2[keepers]

otu_cont2[,1:4][otu_cont2[,1:4] < 0] <- 0

new_cols <- c("as.mock.comm", "jj.mock.comm", "on.mock.comm", "dj.mock.comm")
colnames(otu_cont2) <- new_cols

mock_sams <- subset_samples(pseq, Leaf == 'mock')
mock_sams <- sample_data(mock_sams)

mock_otus <- otu_table(otu_cont2, taxa_are_rows = TRUE)
mocks <- merge_phyloseq(mock_otus, mock_sams, full_tax2)

# sample_sums(mocks)
# set.seed(777)
# mocks_rare <- rarefy_even_depth(mocks, sample.size = 2970, replace = FALSE)

# Okay, run a PERMANOVA and see if they're different

# mock_bray <- phyloseq::distance(mocks_rare, method = 'bray')
# mock_sam_df <- data.frame(sample_data(mocks_rare))
# adonis2(mock_bray ~ Month, data = mock_sam_df)
# This doesn't seem informative with only 3 Df
# Try ordination?

# p_mock <- plot_ordination(mocks_rare, ordinate(mocks_rare, method = 'NMDS', 
#                                                distance = 'bray', trymax = 999), 
#                           'sites', color = 'Month') + stat_ellipse() + 
#   theme_classic()
# p_mock

# I think what I'll need to do is include all of the samples merged by plant
# and add these mocks back in, then ordinate to see if they are indeed 
# different. Also, will need to come up with a more appropriate test
# for significance between them. Look to the literature about what to 
# do with MC data, or on the EMP website

# 2020-03-16
# Will try to observe differences between MCs between runs with 
# stacked barplots

# rank_names(mocks_rare)
# colnames(tax_table(mocks_rare)) =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# 
# mock_rel <- transform_sample_counts(mocks_rare, function(x) x/2970)
# mock_trim <- filter_taxa(mock_rel, function(x) mean(x) > 1e-2, TRUE)
# 
# p3 <- plot_bar(mock_trim, "Family", fill = "Genus")
# p3 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + 
#   ylab("Relative abundance") + facet_wrap(~ Month)
# theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), 
#       panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# colSums(otu_table(mocks_rare))

### Merging the MC pseq object and the full dataset objects
## Will use pseq_merge3 from line 284 above

pseq_merge_mocks <- merge_phyloseq(pseq_merge3, mocks)

sample_data(pseq_merge_mocks)[201:204,2] <- 'mock'

pseq_merge_mocks <- subset_samples(pseq_merge_mocks, Month != 'Apr')

sample_data(pseq_merge_mocks)$Month <- factor(sample_data(pseq_merge_mocks)$Month, levels = 
                                                c('May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov',
                                                  'Dec', 'Jan', 'mock'))

sort(sample_sums(pseq_merge_mocks))

set.seed(777)
pseq_mocks_rare <- rarefy_even_depth(pseq_merge_mocks, sample.size = 2996, 
                                     replace = FALSE)

p3 <- plot_ordination(pseq_mocks_rare, ordinate(pseq_mocks_rare, 'NMDS', distance = 'bray', trymax = 100), 
                      'sites', color = 'Month', title = 'Stress = 0.11') +
  stat_ellipse() + theme_classic()

p3 

ggsave('../R/Figures/ord_supp_v1.pdf', width = d_w, height = d_h, units = 'in')

### 2020/05/08
## Working on getting everything organized and analyzed again 
## Still need to incorporate the within-plant analysis
## I think including the mock communities in ordination is the
  # way to go...

### 2020-05-12
## Import the zotu tables and compare to the OTU tables
  # Going to hold off on the comparison for now
  # Can run correlations between OTUs and ZOTUs of 
    # alpha and beta diversity for each sample
    # but don't need to do it right now

## Processing the within-plant samples:

# Not sure what to do about zotus/clustering. 
  # Maybe just keep separate since they were collected and
  # and sequenced separately

wp_otu <- 'within_p_zotutab2.txt'
wp_map <- 'within_p_map_file.txt'
wp_tax <- 'within_p_zotus3.sintax'

wp_full_otu <- read.delim(wp_otu, header = TRUE, row.names = 1)
wp_full_tax <- read.delim(wp_tax, header = FALSE, row.names = 1)
wp_tax_mat <- as.matrix(wp_full_tax)

wp_full_otu2 <- otu_table(wp_full_otu, taxa_are_rows = TRUE)
wp_full_tax2 <- tax_table(wp_tax_mat)
wp_full_map <- import_qiime_sample_data(wp_map)

wp_pseq <- merge_phyloseq(wp_full_otu2, wp_full_tax2, wp_full_map)

rank_names(wp_pseq)

colnames(tax_table(wp_pseq)) = c('Kingdom', 'Phylum', 'Class',
                                 'Order', 'Family', 'Genus',
                                 'Species')

sort(sample_sums(wp_pseq))

set.seed(777)

wp_rare <- rarefy_even_depth(wp_pseq, sample.size = 5000, replace = FALSE)

# Ordination:

p1 <- plot_ordination(wp_rare, ordinate(wp_rare, 'RDA', distance = 'bray', trymax = 100), 
                      'sites', color = 'Plant_compartment') + stat_ellipse() + theme_classic()

p1

# Rhizome is different, not much else
library('PCDimension') # Use this for broken stick calcs
  # Need to review based on the # of PC dimensions

# Taxonomy plot
wp_ra <- transform_sample_counts(wp_rare, function(x) x/20000)
wp_filt <- filter_taxa(wp_ra, function(x) mean(x) > 2e-3, prune = TRUE)

p2 <- plot_bar(wp_filt, 'Family', fill = 'Genus')
p2 + geom_bar(aes(color = Genus, fill = Genus), stat = 'identity',
              position = 'stack') + ylab('Relative abundance') + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), 
  panel.background = element_blank(), axis.line = 
    element_line(colour = "black"), legend.position = 'none') + 
  facet_wrap(~ Plant_compartment, ncol = 3, nrow = 1)

# Will adjust some more later

tax_table(subset_taxa(wp_rare, Genus %in% 'Derxomyces'))
tax_table(subset_taxa(wp_rare, Genus %in% 'Tremella'))
tax_table(subset_taxa(wp_rare, Genus %in% 'Hymenoscyphus'))

# Alpha diversity plots with hill numbers

wp_otu <- t(otu_table(wp_rare))
hill <- renyi(wp_otu, scales = c(0,1,2), hill = TRUE)
hill2 <- tibble::rownames_to_column(hill, 'SampleID')

library(reshape2)
library(plyr)

hill2 <- melt(hill2, id.vars = 'SampleID')

colnames(hill2) <- c('SampleID', 'Hill_number', 'Hill_value')

hill2['Compartment'] <- as.factor(sub('(\\w+)\\.\\d\\.\\d', '\\1', hill2$SampleID))
hill2['Compartment'] <- mapvalues(hill2$Compartment, from = c("Le", "Ra", "Rh"), 
          to = c("Leaf", "Rachis", "Rhizome"))


p3 <- ggplot(hill2, aes(x=Compartment, y=Hill_value, fill=Compartment)) +
  geom_boxplot() + facet_grid(Hill_number ~., scales = 'free_y')
p3

# Stats
hill_num2 <- subset(hill2, Hill_number == 0)
hill_an2 <- aov(Hill_value ~ Compartment, data = hill_num2)
summary(hill_an2)
plot(TukeyHSD(aov(hill_num2$Hill_value~hill_num2$Compartment)))

# # The base plot
# h <- hillplot <- ggplot(Hillfull2, aes(x=Month, y=Hill.Number, fill=Month)) +
#   geom_boxplot() + facet_grid(q.value ~., scales = "free_y")
# h
# # Change the y-axis labels ?
# h <- hillplot <- ggplot(Hillfull2, aes(x=Month, y=Hill.Number, color=Month)) +
#   geom_boxplot() + facet_grid(q.value ~., scales = "free_y", switch = "y",
#                               labeller = as_labeller(c(Hill_0 = "Richness (OTUs)", Hill_1 = "Exp. Shannon Entropy", 
#                                                        Hill_2 = "Inv. Simpson Index"))) +
#   ylab(NULL) + theme(strip.background = element_blank(), strip.placement = "outside",
#                      strip.text = element_text(size = 12), axis.title = element_text(size = 12), 
#                      axis.text = element_text(size = 12, colour = "black"), legend.position = "none")
# h
# # Saved as Hill2.pdf
# 
# ### Stats on Hill numbers
# 
# # Separate by Hill number
# Hillzero <- dplyr::filter(Hillfull2, q.value == "Hill_0")
# Hillone <- dplyr::filter(Hillfull2, q.value == "Hill_1")
# Hilltwo <- dplyr::filter(Hillfull2, q.value == "Hill_2")
# 
# # Repeated measures anova: https://www.gribblelab.org/stats/notes/RepeatedMeasuresANOVA.pdf
# anzero <- aov(Hill.Number ~ Month + Error(Plant/Month), data = Hillzero)
# summary(anzero)
# 
# anone <- aov(Hill.Number ~ Month + Error(Plant/Month), data = Hillone)
# summary(anone)
# 
# antwo <- aov(Hill.Number ~ Month + Error(Plant/Month), data = Hilltwo)
# summary(antwo)
# 
# plot(TukeyHSD(aov(Hillzero$Hill.Number~Hillzero$Month)), las = 1)











