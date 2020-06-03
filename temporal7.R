## R Code for Younginger, Stewart, Balkan and Ballhorn
## "Stable coexistence or competitive exclusion? 
## Fern endophytes shift our understanding of microbial temporal turnover"

## Outline separate components of the project
## Outline of the separate analyses and plots
## 

setwd('/Users/brett/Documents/Temporal_turnover/Informatics')

# Load packages
library(phyloseq) # 1.32.0
library(vegan) # 2.5-6
# library(lattice) # 0.20-41
# library(permute) # 0.9-5
library(ggplot2) # 3.3.0
library(dplyr) # 0.8.5
library(reshape2) # 1.4.4
library(plyr) # 1.8.6
library(dunn.test) # 1.3.5

# Figure dimensions for journal (ESA)
prop = 5/7      
s_w = 3.00
s_h = s_w*prop
h_w = 4.50
h_h = h_w*prop
d_w = 6.00
d_h = d_w*prop

# Other variables for figure consistency
# outlier size
outlier <- 0.5
# annotation hjust
hjust <- 'left'
# annotation size
size = 2

# Custom theme
theme_w_legend <- function () {
  theme(axis.text = element_text(color = 'black', size = 6),
        axis.line = element_line(color = 'black'),
        axis.title = element_text(size = 8),
        legend.key.height = unit(0.1, 'in'),
        panel.background = element_blank(),
        strip.text = element_text(size = rel(0.5)),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 5),
        legend.margin = margin(l = 0, r = 0, unit = 'pt'),
        legend.key.width = unit(0.1, unit = 'in'), 
        legend.key = element_rect(fill = NA))
}

# Load files

# ZOTUs (ASVs or ESVs):
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

# Separate by month and remove controls (will add back later)
apr_may <- subset_samples(pseq, Month == 'April' | Month == 'May')
jun_jul <- subset_samples(pseq, Month == 'June' | Month == 'July')
aug_sep <- subset_samples(pseq, Month == 'August' | Month == 'September')
oct_nov <- subset_samples(pseq, Month == 'October' | Month == 'November')
dec_jan <- subset_samples(pseq, Month == 'December' | Month == 'January')

# Check sample sizes and make sure the controls are in respective months
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

# All samples accounted for
# Extract the otu tables and sum controls,
# then substract from rest of samples:

# apr_may
apr_may_otu <- otu_table(apr_may)
am_df <- as.data.frame(apr_may_otu)
am_minus <- am_df[,1:159] - am_df[,160]
am_minus[,1:159][am_minus[,1:159] < 0] <- 0
am_clean <- otu_table(am_minus, taxa_are_rows = TRUE)
am_sam <- subset_samples(sample_data(apr_may), Plant != 'control')
am_tax <- tax_table(apr_may)
am_new <- merge_phyloseq(am_sam, am_tax, am_clean)

# jun_jul
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

# aug_sep
aug_sep_df <- as.data.frame(otu_table(aug_sep))
aug_sep_df['all_cont'] <- aug_sep_df['as.mock.comm'] + aug_sep_df['as.pcr.cont']
as_minus <- aug_sep_df[,1:155] - aug_sep_df[,158]
as_minus2 <- as_minus[,1:155]
as_minus2[,1:155][as_minus2[,1:155] < 0] <- 0
as_clean <- otu_table(as_minus2, taxa_are_rows = TRUE)
as_sam <- subset_samples(sample_data(aug_sep), Plant != 'control')
as_tax <- tax_table(aug_sep)
as_new <- merge_phyloseq(as_clean, as_sam, as_tax)

# oct_nov
on_df <- as.data.frame(otu_table(oct_nov))
on_df['all_cont'] <- on_df['on.ex.cont.1'] + on_df['on.ex.cont.2'] + on_df['on.mock.comm']
on_minus <- on_df[,1:153] - on_df[,157]
on_minus2 <- on_minus[,1:153]
on_minus2[,1:153][on_minus2[,1:153] < 0] <- 0
on_clean <- otu_table(on_minus2, taxa_are_rows = TRUE)
on_tax <- tax_table(oct_nov)
on_sam <- subset_samples(sample_data(oct_nov), Plant != 'control')
on_new <- merge_phyloseq(on_clean, on_tax, on_sam)

# dec_jan
dj_df <- as.data.frame(otu_table(dec_jan))
dj_df['all_cont'] <- dj_df[,157] + dj_df[,158] + dj_df[,159] + dj_df[,160]
dj_minus <- dj_df[,1:156] - dj_df[,161]
dj_minus[,1:156][dj_minus[,1:156] < 0] <- 0
dj_clean <- otu_table(dj_minus, taxa_are_rows = TRUE)
dj_tax <- tax_table(dec_jan)
dj_sam <- subset_samples(sample_data(dec_jan), Plant != 'control')
dj_new <- merge_phyloseq(dj_clean, dj_tax, dj_sam)

# Recombine full dataset without controls
pseq_clean <- merge_phyloseq(am_new, jj_new, as_new, on_new, dj_new)

# Declare taxonomic ranks
rank_names(pseq_clean)
colnames(tax_table(pseq_clean)) =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Merge sampling replicates by plant
pseq_merge <- merge_samples(pseq_clean, 'Plant')

# Reformat sample data
sam_mer <- as.data.frame(sample_data(pseq_merge))
sam_mer_rows <- rownames(sam_mer)
sam_mer['New_SampleID'] <- sam_mer_rows
sam_mer['New_Month'] <- sub('\\d+\\.(\\w+)', '\\1', sam_mer$New_SampleID)
sam_mer['New_Plant'] <- as.factor(sub('(\\d+)\\.\\w+', '\\1', sam_mer$New_SampleID))
sam_mer <- sam_mer[,5:7]
new_names <- c('SampleID', 'Month', 'Plant')
colnames(sam_mer) <- new_names

# Add updated sample data to the tax and otu files
mer_tax <- tax_table(pseq_merge)
mer_otu <- otu_table(pseq_merge)
pseq_merge2 <- merge_phyloseq(mer_otu, mer_tax, sam_mer)

## Remove OTUs with < 2 reads
pseq_merge3 <- prune_taxa(taxa_sums(pseq_merge2) > 1, pseq_merge2)

# Rarefy to even sampling depth
set.seed(555)
sort(sample_sums(pseq_merge3))
pseq_rare <- rarefy_even_depth(pseq_merge3, sample.size = 15000, replace = FALSE)
# Lost 11 samples and 45 OTUs

# Remove the seq month April hereafter due to low read numbers, significant loss 
  # of samples after rarefaction. Likely the result of low fungal template 
  # in early-emerging plant foliar tissue

# Rarefaction curves
MAY <- subset_samples(pseq_rare, Month == 'May')
JUN <- subset_samples(pseq_rare, Month == 'Jun')
JUL <- subset_samples(pseq_rare, Month == 'Jul')
AUG <- subset_samples(pseq_rare, Month == 'Aug')
SEP <- subset_samples(pseq_rare, Month == 'Sep')
OCT <- subset_samples(pseq_rare, Month == 'Oct')
NOV <- subset_samples(pseq_rare, Month == 'Nov')
DEC <- subset_samples(pseq_rare, Month == 'Dec')
JAN <- subset_samples(pseq_rare, Month == 'Jan')

MAY <- prune_taxa(taxa_sums(MAY) > 0, MAY)
JUN <- prune_taxa(taxa_sums(JUN) > 0, JUN)
JUL <- prune_taxa(taxa_sums(JUL) > 0, JUL)
AUG <- prune_taxa(taxa_sums(AUG) > 0, AUG)
SEP <- prune_taxa(taxa_sums(SEP) > 0, SEP)
OCT <- prune_taxa(taxa_sums(OCT) > 0, OCT)
NOV <- prune_taxa(taxa_sums(NOV) > 0, NOV)
DEC <- prune_taxa(taxa_sums(DEC) > 0, DEC)
JAN <- prune_taxa(taxa_sums(JAN) > 0, JAN)

May_otu <- as.data.frame(otu_table(MAY))
Jun_otu <- as.data.frame(otu_table(JUN))
Jul_otu <- as.data.frame(otu_table(JUL))
Aug_otu <- as.data.frame(otu_table(AUG))
Sep_otu <- as.data.frame(otu_table(SEP))
Oct_otu <- as.data.frame(otu_table(OCT))
Nov_otu <- as.data.frame(otu_table(NOV))
Dec_otu <- as.data.frame(otu_table(DEC))
Jan_otu <- as.data.frame(otu_table(JAN))

May_otu_sac <- specaccum(May_otu, permutations = 999, method = 'rarefaction')
Jun_otu_sac <- specaccum(Jun_otu, permutations = 999, method = 'rarefaction')
Jul_otu_sac <- specaccum(Jul_otu, permutations = 999, method = 'rarefaction')
Aug_otu_sac <- specaccum(Aug_otu, permutations = 999, method = 'rarefaction')
Sep_otu_sac <- specaccum(Sep_otu, permutations = 999, method = 'rarefaction')
Oct_otu_sac <- specaccum(Oct_otu, permutations = 999, method = 'rarefaction')
Nov_otu_sac <- specaccum(Nov_otu, permutations = 999, method = 'rarefaction')
Dec_otu_sac <- specaccum(Dec_otu, permutations = 999, method = 'rarefaction')
Jan_otu_sac <- specaccum(Jan_otu, permutations = 999, method = 'rarefaction')

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

colnames(May_otu_sac_df) <- colids
colnames(Jun_otu_sac_df) <- colids
colnames(Jul_otu_sac_df) <- colids
colnames(Aug_otu_sac_df) <- colids
colnames(Sep_otu_sac_df) <- colids
colnames(Oct_otu_sac_df) <- colids
colnames(Nov_otu_sac_df) <- colids
colnames(Dec_otu_sac_df) <- colids
colnames(Jan_otu_sac_df) <- colids

May_otu_sac_df$month <- 'May'
Jun_otu_sac_df$month <- 'Jun'
Jul_otu_sac_df$month <- 'Jul'
Aug_otu_sac_df$month <- 'Aug'
Sep_otu_sac_df$month <- 'Sep'
Oct_otu_sac_df$month <- 'Oct'
Nov_otu_sac_df$month <- 'Nov'
Dec_otu_sac_df$month <- 'Dec'
Jan_otu_sac_df$month <- 'Jan'

full_curve <- rbind(May_otu_sac_df, Jun_otu_sac_df, Jul_otu_sac_df, Aug_otu_sac_df, Sep_otu_sac_df, Oct_otu_sac_df, Nov_otu_sac_df, Dec_otu_sac_df, Jan_otu_sac_df)

full_curve$month <- factor(full_curve$month, levels = c('May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'))

p1 <- ggplot(full_curve, aes(x = plant, y = richness, colour = month)) + geom_line() + 
  xlim(1,20) + theme_w_legend() + xlab('Plants') + ylab('Richness')

p1

# ggsave('../R/Figures/rare_v2.pdf', width = s_w, height = s_h, units = 'in')

# Ordination
pseq_rare_noApr <- subset_samples(pseq_rare, Month != 'Apr')

sample_data(pseq_rare_noApr)$Month <- factor(sample_data(pseq_rare_noApr)$Month, 
                                             levels = c('May', 'Jun', 'Jul','Aug', 'Sep', 
                                                        'Oct', 'Nov', 'Dec', 'Jan'))

p2 <- plot_ordination(pseq_rare_noApr, 
                      ordinate(pseq_rare_noApr, 'NMDS', distance = 'bray', trymax = 100), 
                      'sites', color = 'Month') + stat_ellipse() + theme_w_legend() +
  geom_point(size = 0.5)

# Ordination stats
perm1 <- phyloseq::distance(pseq_rare_noApr, method = 'bray')
perm_df <- data.frame(sample_data(pseq_rare_noApr))
adonis2(perm1 ~ Month, strata = Plant, data = perm_df) # Difference between months, constrained by plant
beta <- betadisper(perm1, perm_df$Month)
permutest(beta)
# label1 <- c("PERMANOVA\nPseudo-F", "= 9.61 p < 0.001")

p2 +
  annotate('text', x = -1.5, y = 2.4, hjust = hjust, size = size,
           label = expression(paste("Pseudo-",F['8,170']," = 9.61 p < 0.001"))) + 
  annotate('text', x = -1.5, y = 2.55, hjust = hjust, size = size, label = "PERMANOVA") + 
  annotate('text', x = -1.5, y = 2.05, hjust = hjust, size = size, label = "PERMDISP") + 
  annotate('text', x = -1.5, y = 1.9, hjust = hjust, size = size,
           label = expression(paste("Pseudo-",F['8,170']," = 33.33 p < 0.001"))) + 
  annotate('text', x = -1.5, y = 2.8, hjust = hjust, size = size, 
           fontface = 'bold', label = "Stress = 0.15")

# ggsave('../R/Figures/ord_v3.pdf', width = h_w, height = h_h, units = 'in')

# Processing the mock community seq control samples for inclusion in
  # supplemental plot
controls <- subset_samples(pseq, Plant == 'control')
# nsamples(controls)
# sample_data(controls)
# colnames(otu_table(controls))

# Substract the extraction and pcr controls from each run, leaving 
  # only seq controls
otu_cont <- as.data.frame(otu_table(controls))
otu_cont['as.mock.minus'] <- otu_cont['as.mock.comm'] - otu_cont['as.pcr.cont']
otu_cont['jj.mock.minus'] <- otu_cont['jj.mock.comm'] - otu_cont['jj.ex.cont']
otu_cont['on.ex.sum'] <- otu_cont['on.ex.cont.1'] + otu_cont['on.ex.cont.2']
otu_cont['on.mock.minus'] <- otu_cont['on.mock.comm'] - otu_cont['on.ex.sum']
otu_cont['dj.ex.sum'] <- otu_cont['dj.ex.cont.2'] + otu_cont['dj.ex.cont.3'] + 
  otu_cont['dj.pcr.cont']
otu_cont['dj.mock.minus'] <- otu_cont['dj.mock.comm'] - otu_cont['dj.ex.sum']

# Keep only relevant columns
otu_cont2 <- otu_cont[,13:18]
keepers <- c("as.mock.minus", "jj.mock.minus", "on.mock.minus", "dj.mock.minus")
otu_cont2 <- otu_cont2[keepers]

# Impute zeros for missing data
otu_cont2[,1:4][otu_cont2[,1:4] < 0] <- 0

# Rename columns
new_cols <- c("as.mock.comm", "jj.mock.comm", "on.mock.comm", "dj.mock.comm")
colnames(otu_cont2) <- new_cols

# Subset out only seq control samples and remerge otu, tax and sample data
mock_sams <- subset_samples(pseq, Leaf == 'mock')
mock_sams <- sample_data(mock_sams)
mock_otus <- otu_table(otu_cont2, taxa_are_rows = TRUE)
mocks <- merge_phyloseq(mock_otus, mock_sams, full_tax2)

# Merge the seq controls variable and the full dataset variable
  # excluding extraction and pcr controls
pseq_merge_mocks <- merge_phyloseq(pseq_merge3, mocks)
sample_data(pseq_merge_mocks)[201:204,2] <- 'mock'

# Remove April again
pseq_merge_mocks <- subset_samples(pseq_merge_mocks, Month != 'Apr')
sample_data(pseq_merge_mocks)$Month <- factor(sample_data(pseq_merge_mocks)$Month, levels = 
                                                c('May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov',
                                                  'Dec', 'Jan', 'mock'))

# Rarefy to even depth
sort(sample_sums(pseq_merge_mocks))
set.seed(555)
pseq_mocks_rare <- rarefy_even_depth(pseq_merge_mocks, sample.size = 2996, 
                                     replace = FALSE)
# Lost 225 OTUs

# Supplemental ordination plot with seq controls
p3 <- plot_ordination(pseq_mocks_rare, ordinate(pseq_mocks_rare, 'NMDS', 
                                                distance = 'bray', trymax = 100), 
                      'sites', color = 'Month') + stat_ellipse() + theme_w_legend()

p3 + annotate('text', x = 0, y = -3, hjust = hjust, size = size, fontface = 'bold', 
              label = 'Stress = 0.11')

#ggsave('../R/Figures/ord_supp_v2.pdf', width = h_w, height = h_h, units = 'in')

# Processing the within-plant analysis (separate sampling and analysis)

# Load files
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

# Declare taxonomic ranks
rank_names(wp_pseq)
colnames(tax_table(wp_pseq)) = c('Kingdom', 'Phylum', 'Class',
                                 'Order', 'Family', 'Genus',
                                 'Species')

# Rarefy to even sampling depth
sort(sample_sums(wp_pseq))
set.seed(555)
wp_rare <- rarefy_even_depth(wp_pseq, sample.size = 5000, replace = FALSE)
# Lost 83 OTUs

# Ordination
p4 <- plot_ordination(wp_rare, 
                      ordinate(wp_rare, 'PCoA', distance = 'bray', trymax = 100), 
                      'sites', color = 'Plant_compartment') + 
  stat_ellipse() + theme_w_legend()

p4

# ggsave('../R/Figures/wp_ord_supp_v2.pdf', width = h_w, height = h_h, units = 'in')

# Taxonomy plot
# Transform to relative abundance and filter low abundance taxa
wp_ra <- transform_sample_counts(wp_rare, function(x) x/20000)
wp_filt <- filter_taxa(wp_ra, function(x) mean(x) > 1e-3, prune = TRUE)

p5 <- plot_bar(wp_filt, 'Class', fill = 'Genus')
p5 + geom_bar(aes(color = Genus, fill = Genus), stat = 'identity',
              position = 'stack') + ylab('Relative abundance') + 
  theme_w_legend() + theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                                      vjust = 1)) +
  facet_wrap(~ Plant_compartment, ncol = 3, nrow = 1)

# ggsave('../R/Figures/wp_tax_supp_v2.pdf', width = h_w, height = h_h, units = 'in')

# Alpha diversity plots
# Remove seq month April
pseq_rare_noApr = subset_samples((pseq_rare), Month != 'Apr')
pseq_rare_noApr = prune_taxa(taxa_sums(pseq_rare_noApr) > 0, pseq_rare_noApr)

# Hill numbers
pseq_rare_noApr_otu = otu_table(pseq_rare_noApr)
hill <- renyi(pseq_rare_noApr_otu, scales = c(0,1,2), hill = TRUE)
hill <- tibble::rownames_to_column(hill, 'SampleID')

# Melt dataframe to tidy
hill <- melt(hill, id.vars = 'SampleID')
colnames(hill) = c('SampleID', 'Hill_number', 'Hill_value')

# Adjust grouping variables in dataset
hill$Plant <- as.factor(sub('(\\d+)\\.\\w+', '\\1', hill$SampleID))
hill$Month <- as.factor(sub('\\d+\\.(\\w+)', '\\1', hill$SampleID))
hill$Hill_number <- mapvalues(hill$Hill_number, from = c(0,1,2), 
                              to = c('zero', 'one', 'two'))

hill$Month <- factor(hill$Month, 
                     levels = c('May', 'Jun', 'Jul', 'Aug', 'Sep', 
                                'Oct', 'Nov', 'Dec', 'Jan'))

# Plot all Hill numbers (for supplemental)
p6 <- ggplot(hill, aes(x = Month, y = Hill_value, fill = Month)) + 
  geom_boxplot(outlier.size = outlier) + ylab('Hill value') +
  facet_grid(Hill_number ~., scales = 'free_y', 
             labeller = as_labeller(c('zero'='Richness (OTUs)', 
                                      'one'= 'Exp. Shannon\n Entropy', 
                                      'two'='Inv. Simpson\n Index'))) +
  theme_w_legend() + theme(legend.position = 'none')

p6

# ggsave('../R/Figures/hill_supp_v2.pdf', width = h_w, height = h_h, units = 'in')

# Stats on Hill where q = 0 
hill2 <- subset(hill, Hill_number == 'zero')

hillaov <- aov(Hill_value ~ Month + Error(Plant/Month), data = hill2)
# Get an error: Error() model is singular
# This is likely because Jul only has 19 individuals
  # will drop this individual from the analysis

jul_only <- subset(hill2, Month == 'Jul')
nrow(jul_only)
jul_only # Plant 6 to remove

no_six <- subset(hill2, Plant != 6)
nrow(no_six)

hillaov2 <- aov(Hill_value ~ Month + Error(Plant/Month), data = no_six)
summary(hillaov2)

# Check parametric assumptions
qqnorm(no_six$Hill_value)
qqline(no_six$Hill_value)
hist(no_six$Hill_value)
bartlett.test(Hill_value~Month, data = no_six)

# Log transformation
no_six$log_val <- log(no_six$Hill_value)
qqnorm(no_six$log_val)
qqline(no_six$log_val)
bartlett.test(log_val~Month, data = no_six)
hist(no_six$log_val)

# Run repeated-measures ANOVA with log-transformed data
hillaov3 <- aov(log_val ~ Month + Error(Plant/Month), data = no_six)
summary(hillaov3)

Tuk <- TukeyHSD(aov(log_val ~ Month, data = no_six), las = 1)

# Compare with Kruskal
kruskal.test(Hill_value ~ Month, data = no_six)
# Same result

## Capture non-significant post-hoc values for plot annotation
Tuk2 <- Tuk$Month
Tuk2 <- as.data.frame(Tuk2)
Tuk2$p_adj <- Tuk2$`p adj`
non_sig <- subset(Tuk2, p_adj >= 5e-2)
non_sig

# Stats annotations
labels <- c(' A', ' A', ' B', ' C', ' C', ' D', 'D,B', 'D,B', ' D')
x_locs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
tapply(hill2$Hill_value, hill2$Month, summary)  # To figure out the value of the upper quartile
y_locs <- c(92.25, 91.50, 11.50, 7.00, 7.75, 27.50, 15.75, 20.75, 24.25)

# Alpha diversity plot
p7 <- ggplot(hill2, aes(x = Month, y = Hill_value, fill = Month)) + 
  geom_boxplot(outlier.size = outlier) + 
  theme_w_legend() + theme(legend.position = 'none') + ylab('Richness (OTUs)') +
  annotate('text', x = 3, y = 93, hjust = hjust, size = size, 
           label = 'Repeated measures ANOVA') + 
  annotate('text', x = 3, y = 85, hjust = hjust, size = size, 
           label = expression(paste(F["8,144"]," = 103.5, p < 0.001"))) + 
  annotate('text', x = x_locs, y = y_locs, label = labels, vjust = -0.5, 
           hjust = -0.1, size = size)

p7

# ggsave('../R/Figures/hill_v5.pdf', width = s_w, height = s_h, units = 'in')

# Taxonomy plot
# Convert to relative abundance
pseq_rel <- transform_sample_counts(pseq_rare_noApr, function(x) x/3e+05)

# Supplemental plot without filtering
sample_data(pseq_rel)$Month <- factor(sample_data(pseq_rel)$Month, 
                                      levels = c('May', 'Jun', 'Jul',
                                                 'Aug', 'Sep', 'Oct',
                                                 'Nov', 'Dec', 'Jan'))

p8 <- plot_bar(pseq_rel, x = 'Class', fill = 'Genus')

p8 + geom_bar(aes(fill = Genus, color = Genus), stat = 'identity', 
              position = 'stack') + facet_wrap(~Month, nrow = 3) + 
  theme_w_legend() + 
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  ylab('Relative abundance')

# ggsave('../R/Figures/tax_supp_V2.pdf', width = d_w, height = d_h, units = 'in')

# Final plot with filtering
# Filter low-abundance taxa
pseq_rel_filt <- filter_taxa(pseq_rel, function(x) mean(x) > 2e-5, prune = TRUE)

p9 <- plot_bar(pseq_rel_filt, x = 'Class', fill = 'Genus') + 
  geom_bar(aes(fill = Genus, color = Genus), stat = 'identity', 
           position = 'stack') + facet_wrap(~Month, nrow = 3) + 
  theme_w_legend() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1)) + 
  ylab('Relative abundance')

p9

# ggsave('../R/Figures/tax_V2.pdf', width = h_w, height = h_h, units = 'in')

# Gametophyte bioassay results
# Import data
gamete <- read.csv('../Gametophytes/Nov_time1.csv', header = T, sep = "\t")

# Adjust data
gamete$Sample_ID2 <- sub('(\\w+_\\d+)_\\d', '\\1', gamete$SampleID)

# Add taxonomic names to accession numbers
gamete$SampleID[gamete$Sample_ID2 == 'BSY_163'] <- 'Flagellospora'
gamete$SampleID[gamete$Sample_ID2 == 'BSY_137'] <- 'Plectania'
gamete$SampleID[gamete$Sample_ID2 == 'Cont_1'] <- 'Control'
gamete$SampleID[gamete$Sample_ID2 == 'Cont_2'] <- 'Control'

# Stats on affects of fungal inocula on gametophyte growth
week4 <- subset(gamete, Time_elapsed == '4 weeks')
week8 <- subset(gamete, Time_elapsed == '8 weeks')

# Not likely to meet parametric assumptions. Will use a 
  # Kruskal-Wallis rank sum with a Hochberg correction
  # Analyze each timepoint separately, comparing samples within
  # a timepoint, rather than as repeated measures

week4_mod1 <- kruskal.test(Nov_Time1 ~ SampleID, data = week4)
week4_mod1

week8_mod1 <- kruskal.test(Nov_Time1 ~ SampleID, data = week8)
week8_mod1

week8$SampleID <- as.factor(week8$SampleID)
dunn.test(week8$Nov_Time1, week8$SampleID, method = 'hochberg')

f_labels <- data.frame(x = c(1, 1), y = c(-0.25, -0.25), 
                       Time_elapsed = c('4 weeks', '8 weeks'), 
                       label = c('chi^2(2) = 4.00, p = 0.136', 
                                 'chi^2(2) = 17.11, p < 0.001'))

# Impute dummy data into the dataframe to obtain separate annotations
  # in facets with different statistical results
gamete2 <- gamete
gamete2[,5:7] <- NA
newcols <- c("SampleID", "Nov_Time1", "Time_elapsed", "Sample_ID2", "X", "Y", "label")
colnames(gamete2) = newcols
gamete2[1,5] <- 1.50
gamete2[31,5] <- 1.50
gamete2[1,6] <- -0.25
gamete2[31,6] <- -0.25
gamete2[1,7] <- '4.00, p = 0.136'
gamete2[31,7] <- '17.11, p < 0.001'

# Figure out the values of the upper quartile to plot the post-hoc values
tapply(week8$Nov_Time1, week8$SampleID, summary)
gamete2[32,5] <- 1.1
gamete2[33,5] <- 2.1
gamete2[34,5] <- 3.1
gamete2[32,6] <- 0.28
gamete2[33,6] <- 0.34
gamete2[34,6] <- 0.04
gamete2[2,7] <- ' '
gamete2[3,7] <- ' '
gamete2[4,7] <- ' '
gamete2[32,7] <- 'A'
gamete2[33,7] <- 'A'
gamete2[34,7] <- 'B'

# The plot
p10 <- ggplot(gamete2, aes(x=SampleID, y=Nov_Time1, fill = SampleID)) + 
  geom_boxplot(outlier.size = outlier) + facet_grid(.~Time_elapsed) +
  theme_w_legend() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        legend.position = 'none') +
  geom_hline(yintercept = 0, alpha = 0.5) + xlab('Treatment group') + 
  ylab(expression(paste(Delta,' Surface area ', (cm)^2))) + 
  annotate('text', x = 0.75, y = -0.25, 
           label = expression(paste(chi^2,'(2) = ')), 
           hjust = hjust, size = size) + 
  annotate('text', x = 0.725, y = -0.18, label = 'Kruskal-Wallis', 
           hjust = hjust, size = size) + 
  geom_text(data = gamete2, aes(x = X, y = Y, label = label), 
            hjust = hjust, size = size)

p10

# ggsave('../R/Figures/gamete_V4.pdf', width = s_w, height = s_h, units = 'in')

# Competition assay plot
# Load data
comp1 <- read.csv('../Competition_assays/compiled_images.csv', header = TRUE)

# EDA of the following:
# - length:towards
# - area:towards
# - area:away

tow_len <- dplyr::filter(comp1, Direction == 'towards')
awa_are <- dplyr::filter(comp1, Direction == 'away')

ggplot(tow_len, aes(x = Competitor, y = Length)) + geom_boxplot()
ggplot(tow_len, aes(x = Competitor, y = Area)) + geom_boxplot()
ggplot(awa_are, aes(x = Competitor, y = Area)) + geom_boxplot()

# Plot the length-towards (amount of growth towards C. polystichicola competitor):

p11 <- ggplot(tow_len, aes(x = Competitor, y = Length, fill = Competitor)) + 
  geom_boxplot(outlier.size = outlier) + theme_w_legend() +
  scale_x_discrete(breaks = c("C","D","E","F","G","I","J","K","L","M","N"), 
                   labels = c("A","B","C","D","E","F","G","H","I","J","K")) +
  scale_fill_discrete(labels = c("A: Arthrinium arundinis","B: Pezicula sp.",
                                 "C: Xylaria sp.","D: Diaporthe sp.","E: Colletotrichum sp.1",
                                 "F: Hypoxylon rubiginosum","G: Colletotrichum sp.2","H: Phoma sp.",
                                 "I: Plectania milleri","J: Fontanospora sp.","K: Pezicula sp.")) +
  scale_y_continuous(name = "Difference in growth relative\nto C. polystichicola (cm)")

p11

# ggsave('../R/Figures/comp_V2.pdf', width = s_w, height = s_h, units = 'in')
