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
#1. Load files

otu <- 'merged_otutab_V4.txt'
map <- '/Users/brett/Documents/Temporal_turnover/Informatics/map_file16.txt'
tax <- 'merged_otus_V3.sintax'

full_otu <- read.delim(otu, header = TRUE, row.names = 1)
full_tax <- read.delim(tax, header = FALSE, row.names = 1)
tax_mat <- as.matrix(full_tax)

full_otu2 <- otu_table(full_otu, taxa_are_rows = TRUE)
full_tax2 <- tax_table(tax_mat)
full_map <- import_qiime_sample_data(map)

pseq <- merge_phyloseq(full_otu2, full_tax2, full_map)

apr_may <- subset_samples(pseq, Month == 'April' | Month == 'May') # '|' is 'or'

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

sample_sums(pseq_clean)

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

sam_sums <- sample_sums(pseq_merge2)
sort(sam_sums)

## Remove some spurious OTUs, potential contaminants? 

pseq_merge3 <- prune_taxa(taxa_sums(pseq_merge2) > 1, pseq_merge2) # Maybe adjust this later? 
pseq_merge3 <- subset_taxa(pseq_merge3, Genus != 'Candida')

## Good for now. Will come back and clean more after visualizing some 
# taxonomic plots
set.seed(777)
pseq_rare <- rarefy_even_depth(pseq_merge3, sample.size = 15000, replace = FALSE)
# Lost 11 samples and 134 OTUs

# Plot some rarefaction curves

APR <- subset_samples(pseq_rare, Month == 'Apr')
MAY <- subset_samples(pseq_rare, Month == 'May')
JUN <- subset_samples(pseq_rare, Month == 'Jun')
JUL <- subset_samples(pseq_rare, Month == 'Jul')
AUG <- subset_samples(pseq_rare, Month == 'Aug')
SEP <- subset_samples(pseq_rare, Month == 'Sep')
OCT <- subset_samples(pseq_rare, Month == 'Oct')
NOV <- subset_samples(pseq_rare, Month == 'Nov')
DEC <- subset_samples(pseq_rare, Month == 'Dec')
JAN <- subset_samples(pseq_rare, Month == 'Jan')

APR <- prune_taxa(taxa_sums(APR) > 0, APR)
MAY <- prune_taxa(taxa_sums(MAY) > 0, MAY)
JUN <- prune_taxa(taxa_sums(JUN) > 0, JUN)
JUL <- prune_taxa(taxa_sums(JUL) > 0, JUL)
AUG <- prune_taxa(taxa_sums(AUG) > 0, AUG)
SEP <- prune_taxa(taxa_sums(SEP) > 0, SEP)
OCT <- prune_taxa(taxa_sums(OCT) > 0, OCT)
NOV <- prune_taxa(taxa_sums(NOV) > 0, NOV)
DEC <- prune_taxa(taxa_sums(DEC) > 0, DEC)
JAN <- prune_taxa(taxa_sums(JAN) > 0, JAN)

Apr_otu <- as.data.frame(otu_table(APR))
May_otu <- as.data.frame(otu_table(MAY))
Jun_otu <- as.data.frame(otu_table(JUN))
Jul_otu <- as.data.frame(otu_table(JUL))
Aug_otu <- as.data.frame(otu_table(AUG))
Sep_otu <- as.data.frame(otu_table(SEP))
Oct_otu <- as.data.frame(otu_table(OCT))
Nov_otu <- as.data.frame(otu_table(NOV))
Dec_otu <- as.data.frame(otu_table(DEC))
Jan_otu <- as.data.frame(otu_table(JAN))

Apr_otu_sac <- specaccum(Apr_otu, permutations = 999, method = 'rarefaction')
May_otu_sac <- specaccum(May_otu, permutations = 999, method = 'rarefaction')
Jun_otu_sac <- specaccum(Jun_otu, permutations = 999, method = 'rarefaction')
Jul_otu_sac <- specaccum(Jul_otu, permutations = 999, method = 'rarefaction')
Aug_otu_sac <- specaccum(Aug_otu, permutations = 999, method = 'rarefaction')
Sep_otu_sac <- specaccum(Sep_otu, permutations = 999, method = 'rarefaction')
Oct_otu_sac <- specaccum(Oct_otu, permutations = 999, method = 'rarefaction')
Nov_otu_sac <- specaccum(Nov_otu, permutations = 999, method = 'rarefaction')
Dec_otu_sac <- specaccum(Dec_otu, permutations = 999, method = 'rarefaction')
Jan_otu_sac <- specaccum(Jan_otu, permutations = 999, method = 'rarefaction')

Apr_otu_sac_df <- data.frame(Apr_otu_sac$sites, Apr_otu_sac$richness, Apr_otu_sac$sd)
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

colnames(Apr_otu_sac_df) <- colids
colnames(May_otu_sac_df) <- colids
colnames(Jun_otu_sac_df) <- colids
colnames(Jul_otu_sac_df) <- colids
colnames(Aug_otu_sac_df) <- colids
colnames(Sep_otu_sac_df) <- colids
colnames(Oct_otu_sac_df) <- colids
colnames(Nov_otu_sac_df) <- colids
colnames(Dec_otu_sac_df) <- colids
colnames(Jan_otu_sac_df) <- colids

Apr_otu_sac_df$month <- 'Apr'
May_otu_sac_df$month <- 'May'
Jun_otu_sac_df$month <- 'Jun'
Jul_otu_sac_df$month <- 'Jul'
Aug_otu_sac_df$month <- 'Aug'
Sep_otu_sac_df$month <- 'Sep'
Oct_otu_sac_df$month <- 'Oct'
Nov_otu_sac_df$month <- 'Nov'
Dec_otu_sac_df$month <- 'Dec'
Jan_otu_sac_df$month <- 'Jan'

full_curve <- rbind(Apr_otu_sac_df, May_otu_sac_df, Jun_otu_sac_df, Jul_otu_sac_df, Aug_otu_sac_df, Sep_otu_sac_df, Oct_otu_sac_df, Nov_otu_sac_df, Dec_otu_sac_df, Jan_otu_sac_df)

full_curve$month <- factor(full_curve$month, levels = c('Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'))

p1 <- ggplot(full_curve, aes(x = plant, y = richness, colour = month)) + geom_line() + xlim(1,20)
p1 + theme_classic() + xlab('Plants') + ylab('Richness') + scale_fill_discrete(name = 'Months')

# Think instead I should set the xlim to 1 # Will work on this further at a later point
# Also, April is really the outlier due to low sample sizes.
  # I think I should be removing it from the analysis, but will think about it some more

# Ordination

p2 <- plot_ordination(pseq_rare, ordinate(pseq_rare, 'NMDS', distance = 'bray', trymax = 100), 'sites', color = 'Month') +
  stat_ellipse() + theme_classic()

p2
  # Need to organize the months as factors

## Taxonomy

sample_data(pseq_rare)$Month <- factor(sample_data(pseq_rare)$Month, levels = c('Apr', 'May', 'Jun', 'Jul',
                                                                                'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'))

# rela = transform_sample_counts(plant3_rare, function(x) x/11050)
# 
# no_low = filter_taxa(rela, function(x) mean(x) > 1e-3, TRUE)

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

sample_sums(mocks)
set.seed(777)
mocks_rare <- rarefy_even_depth(mocks, sample.size = 2970, replace = FALSE)

# Okay, run a PERMANOVA and see if they're different

mock_bray <- phyloseq::distance(mocks_rare, method = 'bray')
mock_sam_df <- data.frame(sample_data(mocks_rare))
adonis2(mock_bray ~ Month, data = mock_sam_df)
# This doesn't seem informative with only 3 Df
# Try ordination?

p_mock <- plot_ordination(mocks_rare, ordinate(mocks_rare, method = 'NMDS', 
                                               distance = 'bray', trymax = 999), 
                          'sites', color = 'Month') + stat_ellipse() + 
  theme_classic()
p_mock

# I think what I'll need to do is include all of the samples merged by plant
# and add these mocks back in, then ordinate to see if they are indeed 
# different. Also, will need to come up with a more appropriate test
# for significance between them. Look to the literature about what to 
# do with MC data, or on the EMP website

# 2020-03-16
# Will try to observe differences between MCs between runs with 
# stacked barplots

rank_names(mocks_rare)
colnames(tax_table(mocks_rare)) =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

mock_rel <- transform_sample_counts(mocks_rare, function(x) x/2970)
mock_trim <- filter_taxa(mock_rel, function(x) mean(x) > 1e-2, TRUE)

p3 <- plot_bar(mock_trim, "Family", fill = "Genus")
p3 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + 
  ylab("Relative abundance") + facet_wrap(~ Month)
theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), 
      panel.background = element_blank(), axis.line = element_line(colour = "black"))

colSums(otu_table(mocks_rare))

### Merging the MC pseq object and the full dataset objects
## Will use pseq_merge3 from line 284 above

pseq_merge_mocks <- merge_phyloseq(pseq_merge3, mocks)

sample_data(pseq_merge_mocks)[201:204,2] <- 'mock'

set.seed(777)
pseq_mocks_rare <- rarefy_even_depth(pseq_merge_mocks, sample.size = 2970, 
                                     replace = FALSE)

p2 <- plot_ordination(pseq_mocks_rare, ordinate(pseq_mocks_rare, 'NMDS', distance = 'bray', trymax = 100), 'sites', color = 'Month') +
  stat_ellipse() + theme_classic()

p2

### 2020/05/08
## Working on getting everything organized and analyzed again 
## Still need to incorporate the within-plant analysis
## I think including the mock communities in ordination is the
  # way to go...







