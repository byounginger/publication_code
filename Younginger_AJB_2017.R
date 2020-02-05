# Younginger & Ballhorn—American Journal of Botany 2017—Appendix S5

# Contents
# 1 Import and processing
# 2 Rarefaction curves
# 3 Rarefy to even sequencing depth
# 4 Alpha diversity with Hill numbers (Some code adapted from Bálint et al. 2015)
# 5 PERMANOVA/beta diversity
# 6 NMDS
# 7 Taxonomy plots

##########

# 1 Import and processing 

library(ggplot2) # Wickham & Chang, 2016 (ver 2.2.1)
library(phyloseq) # McMurdie & Holmes, 2016 (ver 1.19.1)
library(vegan) # Oksanen et al., 2017  (ver 2.4-2)
library(grid) # Murrell, 2016 (ver 3.3.2)
library(gridExtra) # Auguie & Antonov, 2016 (ver 2.2.1)

# Import json-formatted otu table and associated sample data file
  # This otu table has the number of reads found in the negative control 
  # for each taxon subtracted from each sample, as per Nguyen et al., 2015
biom = "Appendix_S2.biom"
mapping_file = "Appendix_S3.txt"
data = import_biom(biom)
envir = import_qiime_sample_data(mapping_file)
pseq = merge_phyloseq(data,envir)
rank_names(pseq)
colnames(tax_table(pseq)) =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Remove negative control
pseq = subset_samples(pseq, X.SampleID !='control')

# Merge samples by plant
plant = merge_samples(pseq, 'Plant')

# Export the sample data file
write.table(sample_data(plant), file = 'Appendix_S4.txt', sep = '\t')
# Adjust the sample data file following the merge and re-import
pl_tb = 'Appendix_S4.txt'
pl_tb = import_qiime_sample_data(pl_tb)

# Generate otu and tax table for new sample data 
otus = otu_table(plant)
tax = tax_table(plant)

# Remerge
plant2 = merge_phyloseq(otus,tax,pl_tb)

############

# 2 Rarefaction curves

# Split apart by month 
apr_only = subset_samples(plant2, Month == 'Newly emerged')
may_only = subset_samples(plant2, Month == 'One month old')

# Check the number of OTUs per month
apr_non <- prune_taxa(taxa_sums(apr_only) > 0, apr_only)
may_non <- prune_taxa(taxa_sums(may_only) > 0, may_only)
ntaxa(apr_non)
ntaxa(may_non)

# Convert to OTU table
apr_otu = otu_table(apr_only)
may_otu = otu_table(may_only)

# Rarefaction curves
SACmay <- specaccum(may_otu, permutations = 999, conditioned = F)
SACapr <- specaccum(apr_otu, permutations = 999, conditioned = F)
plot(SACmay, ci.type="poly", col="darkgreen", lwd=3, ci.lty=0, ci.col= rgb(0,0,0,0.5), xlab = "Plants", ylab = "OTUs")
plot(SACapr, add = T, ci.type="poly", col="blue", lwd=3, ci.lty=0, ci.col= rgb (.5,.5,.5,0.5))
legend("topleft", inset = c(0,-0.05), y.intersp = 0.2, x.intersp = 0.1, legend = c("Newly emerged", "One month old"), col=c("blue", "dark green"),lty = 1, ncol = 1, cex=1, lwd=2, bty = "n")
# Saved as Figure_1.pdf

##########

# 3 Rarefy to even sequencing depth

# Rarefy
set.seed(711)
sample_sums(plant2)
pt_rare = rarefy_even_depth(plant2, sample.size = 1889)
# Will use these data (pt_rare) throughout the remainder of the analyses

##########

# 4 Alpha diversity with Hill numbers

# Convert rarefied Phyloseq object (pt_rare) to otu tables
apr_rare = subset_samples(pt_rare, Month == 'Newly emerged')
may_rare = subset_samples(pt_rare, Month == 'One month old')
apr_rare_otu = otu_table(apr_rare)
may_rare_otu = otu_table(may_rare)

# Hill numbers
Hillapr = renyi(apr_rare_otu, scales = c(0,1,2), hill = T)
hill.1apr = Hillapr$"0"
hill.2apr = Hillapr$"1"
hill.3apr = Hillapr$"2"
Hillmay = renyi(may_rare_otu, scales = c(0,1,2), hill = T)
hill.1may = Hillmay$"0"
hill.2may = Hillmay$"1"
hill.3may = Hillmay$"2"

# Plots 
months = c('Newly emerged','One month old')
par(mfrow = c(1,3))
boxplot(hill.1apr,hill.1may, notch=F, outline = F,
        xlim=c(0.5,2.5), ylab=NA, main="Hill: q = 0", xaxt="n")
axis(side=1, labels = months, at=c(1,2), lwd=0, las=1, line = -0.5, cex.axis=1.2)
boxplot(hill.2apr, hill.2may, notch=F, outline = F,
        xlim=c(0.5,2.5), ylab=NA, main="Hill: q = 1", xaxt="n")
axis(side=1, labels = months, at=c(1,2), lwd=0, las=1, line = -0.5, cex.axis=1.2)
boxplot(hill.3apr, hill.3may, notch=F, outline = F,
        xlim=c(0.5,2.5), ylab=NA, main="Hill: q = 2", xaxt='n')
axis(side=1, labels = months, at=c(1,2), lwd=0, las=1, line = -0.5, cex.axis=1.2)
# Saved as Figure_2.pdf

##########

# 5 PERMANOVA/beta diversity

# PERMANOVA w/ Bray-Curtis
bray <- phyloseq::distance(pt_rare, method = 'bray')
sampledf <- data.frame(sample_data(pt_rare))
adonis(bray ~ Month, data = sampledf, method = 'bray')

# Group dispersions
beta_bray <- betadisper(bray, sampledf$Month)
permutest(beta_bray)

# Note: the remaining code from section 5 is not included
  # in the final analysis

# PERMANOVA with Morisita-Horn 
horn <- phyloseq::distance(pt_rare, method = 'horn')
adonis(horn ~ Month, data = sampledf, method = 'horn')

# Group dispersions
beta_horn <- betadisper(horn, sampledf$Month)
permutest(beta_horn)

# Filter out low abundance taxa (< 0.5% of total dataset)
rela = transform_sample_counts(pt_rare, function(x) x / sum(x))
no_low = filter_taxa(rela, function(x) mean (x) > 5e-3, TRUE) 

# Convert to OTU tables
apr_no_low = subset_samples(no_low, Month == 'Newly emerged')
may_no_low = subset_samples(no_low, Month == 'One month old')
apr_no_low_otu = otu_table(apr_no_low)
may_no_low_otu = otu_table(may_no_low)

# Rarefaction curves to check for asymptotes
par(mfrow = c(1,1))
SACmay <- specaccum(may_no_low_otu, permutations = 999, conditioned = F)
SACapr <- specaccum(apr_no_low_otu, permutations = 999, conditioned = F)
plot(SACapr, ci.type="poly", col="blue", lwd=3, ci.lty=0, ci.col= rgb (.5,.5,.5,0.5))
plot(SACmay, add = T, ci.type="poly", col="darkgreen", lwd=3, ci.lty=0, ci.col= rgb(0,0,0,0.5), xlab = "Plants", ylab = "OTUs")
legend("topleft", inset = c(0,-0.05), y.intersp = 0.2, x.intersp = 0.1, legend = c("Newly emerged", "One month old"), col=c("blue", "dark green"),lty = 1, ncol = 1, cex=1, lwd=2, bty = "n")

# PERMANOVA w/ Bray-Curtis
bray_no_low <- phyloseq::distance(no_low, method = 'bray')
sampledf_no_low <- data.frame(sample_data(no_low))
adonis(bray_no_low ~ Month, data = sampledf_no_low)

# Group dispersions
beta_bray_no_low <- betadisper(bray_no_low, sampledf_no_low$Month)
permutest(beta_bray_no_low)

# PERMANOVA with Morisita-Horn 
horn_no_low <- phyloseq::distance(no_low, method = 'horn')
adonis(horn_no_low ~ Month, data = sampledf_no_low, method = 'horn')

# Group dispersions
beta_horn_no_low <- betadisper(horn_no_low, sampledf_no_low$Month)
permutest(beta_horn_no_low)

##########

# 6 NMDS

# Transform sample counts to relative abundance
rela = transform_sample_counts(pt_rare, function(x) x / sum (x))

# Ordination
q <- plot_ordination(rela, ordinate(rela, "NMDS", distance = "bray"), "sites", color = "Month") + 
  stat_ellipse() + geom_text(aes(label = X.SampleID), show.legend = FALSE, size = 4, vjust = 1.5) + theme_bw()
q + scale_color_manual(values = c("blue","darkgreen")) + theme(legend.text=element_text(size=12), 
  legend.title=element_text(size=0), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate('text', x=0.9, y=-1.86, label = 'PERMANOVA', hjust = 'left') +
  annotate('text', x=0.9, y=-1.98, parse = T, label = 'R ^ 2 == 0.375', hjust = 'left') +
  annotate('text', x=0.9, y=-2.1, label = 'P < 0.001', hjust = 'left')
# Note: Stress = 0.16
# Saved as Figure_3.pdf

###########

# 7 Taxonomy plots

# Transform to relative abundance within month
rela2 = transform_sample_counts(pt_rare, function(x) x / 37780) # Note: this denominator was determined by summing the read numbers from each sample within each month

# Filter low abundance taxa (< 0.03% of total dataset)
no_low2 = filter_taxa(rela2, function(x) mean (x) > 3e-4, TRUE)

# Family plot
p2 <- plot_bar(no_low2, "Family", fill = "Genus", facet_grid = ~Month)
family_plot <- p2 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + 
  ylab("Relative abundance") + scale_y_continuous(limits = c(0,0.7), expand = c(0, 0)) +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  geom_text(data = data.frame(x = 1.2, y = 0.65, label = c("B",""), Month=c('Newly emerged', 'One month old')), size = 6, aes(x,y,label=label), inherit.aes = F)
family_plot

# Transform sample counts to one
no_low3 = transform_sample_counts(no_low2, function(x) x *20)

# Plant plot
p3 = plot_bar(no_low3, "X.SampleID", fill = "Genus", facet_grid = ~Month)
plant_plot <- p3 + geom_bar(aes(color=Genus, fill=Genus), stat = "identity", position = "stack", width = 0.5) + 
  ylab("Relative abundance") + xlab("Plant") + scale_x_continuous(breaks = 1:20) + scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = -90, vjust = -0.18), axis.ticks.x=element_blank()) +
  geom_text(data = data.frame(x = 0, y = 0.95, nudge_x = -3, label = c("A",""), Month=c('Newly emerged', 'One month old')), size = 6, aes(x,y,label=label), inherit.aes = F)
plant_plot

# Plot combined figures with a shared legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  invisible(combined)
}

full_plot <- grid_arrange_shared_legend(plant_plot, family_plot, ncol = 2, nrow = 1)
# Saved as Figure_4.pdf
