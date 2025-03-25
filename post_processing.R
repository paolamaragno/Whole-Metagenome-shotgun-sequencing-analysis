setwd('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/immubac/metaphlan_humann')

library(vegan)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(psych)
library(plotrix)
library(rstatix)
library(readr)
library(tibble)
library(stringr)
library(RColorBrewer)

# in merged_abundance_table the taxon's relative abundance in %. Since typical shotgun-sequencing-based taxonomic profile is relative 
# (i.e. it does not provide absolute cell counts), clades are hierarchically summed. Each taxonomic level will sum to 100%. That is, the sum of 
# all kingdom-level clades is 100%, the sum of all phylum-level clades is 100%, and so forth.

abundance_table <- read.table('merged_abundance_table.txt', header=T)
rownames(abundance_table) <- abundance_table[,1]
abundance_table[,1] <- NULL
colSums(abundance_table)

metadata <- read.table('../metadata.txt')
colnames(abundance_table) <- rownames(metadata)

abundance_table <- as.data.frame(t(abundance_table))
taxa <- colnames(abundance_table)
meta <- colnames(metadata)

df <- merge(metadata, abundance_table, by = "row.names")

row.names(df) <- df$Row.names
df$Row.names <- NULL

####species alpha####

# Shannon diversity index
shannon_div <- diversity(select(df, all_of(taxa)), index = "shannon")

# Simpson diversity index
simpson_div <- diversity(select(df, all_of(taxa)), index = "simpson")

# Richness (total number of taxa present in each sample)
richness <- rowSums(select(df, all_of(taxa)) > 0)

# Create a data frame to combine the diversity indices
alpha_diversity_df <- data.frame(
  Shannon = shannon_div,
  Simpson = simpson_div,
  Richness = richness
)

alpha_diversity_df$Cluster <- df$Cluster
alpha_diversity_df$Cluster <- factor(alpha_diversity_df$Cluster, levels = c('C1','C2','UC'))
alpha_diversity_df$sample <- rownames(alpha_diversity_df)

## welsh t test on observed richness
alpha_diversity_df_observed <- alpha_diversity_df[,c('Richness', 'Cluster')]

alpha_diversity_df_observed %>% group_by('Cluster') %>% t_test(Richness ~ Cluster, var.equal=FALSE)
# C1-C2, pvalue = 0.023 
##

write.table(alpha_diversity_df, file = '../results_r/alpha_indexes.txt', quote=F, sep='\t',row.names = T, col.names = T)

alpha_diversity_df2 <- alpha_diversity_df %>% gather(key = 'index', value ='value', c(-sample, -Cluster))

alpha_diversity_df2 <- group_by(alpha_diversity_df2, index, Cluster) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = std.error(value, na.rm = TRUE)
  )

alpha_diversity_df2_observed <- alpha_diversity_df2[alpha_diversity_df2$index == 'Richness',]

pdf('../results_r/alpha_diversity.pdf', width = 4, height = 3)
ggplot(alpha_diversity_df2_observed,aes(x = Cluster, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, 
                position=position_dodge(0.05)) +
  geom_point(aes(fill = Cluster), size=4, pch = 21, color='black') +
  theme_classic() +
  theme(aspect.ratio = 2.5, axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#f9c74f", "#f8961e","#f94144")) +
  ylab('Observed n. of Species') 
#stat_pvalue_manual(stat.test,label = "p.adj", tip.length=0.01 )
dev.off()

####species beta####
# Beta diversity analysis: Bray-Curtis distance calculation
beta_dist <- vegdist(select(df, all_of(taxa)), method = "bray")
variance_exp <- ape::pcoa(beta_dist)$values$Relative_eig

beta_dist_table <- as.matrix(beta_dist)
write.table(beta_dist_table, file = '../results_r/bray_curtis.txt', quote=F, sep='\t',row.names = T, col.names = T)

# PCoA analysis
pcoa_results <- cmdscale(beta_dist, eig = TRUE, k = 2)  # Perform PCoA with 2 components
pcoa_df <- data.frame(
  MDS1 = pcoa_results$points[, 1],
  MDS2 = pcoa_results$points[, 2],
  cluster = df$Cluster  # Add study_condition for coloring
)

# Set the order of StudyCondition factor
pcoa_df$cluster <- factor(pcoa_df$cluster, levels = c("C1", "C2", "UC"))

# Define color palette
colors <- c("C1" = "#f9c74f", "C2" = "#f8961e", "UC" = "#f94144")

# Perform PERMANOVA (adonis2) to compare groups based on Bray-Curtis distance
set.seed(1004)
adonis_results <- adonis2(beta_dist ~ df$Cluster, permutations = 999)
# 0.001

# PCoA plot with jittered points and p-value
pdf('../results_r/PCA_based_on_taxa.pdf', width = 6, height = 5)
ggplot(pcoa_df, aes(x = MDS1, y = MDS2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +  # Add points with opacity
  scale_color_manual(values = colors) +
  #stat_ellipse(aes(group = StudyCondition), level = 0.95, type = "t", size = 1) +  # Add ellipses
  labs(title = "PCoA of Beta Diversity (Bray-Curtis)", x = "MDS1", y = "MDS2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  annotate("text", x = Inf, y = Inf, label = paste("p = ", round(adonis_results$`Pr(>F)`[1], 3)),
           hjust = 2, vjust = 1.1, size = 5, color = "black")  # Add p-value from PERMANOVA
dev.off()

fit_p <- hclust(beta_dist, method="average")
pdf('../results_r/hierarchical_tree.pdf', width = 5, height = 5)
plot(fit_p, hang = -1, cex = 1, main = "Bray-Curtis Dissimilarity", 
     ylab = "Bray-Curtis Dissimilarity", xlab = "ID")
dev.off()


####barplot relative abundances####
# Extract the abundance data (only taxa columns) 
abundance_data <- df[, taxa]

## keep only species
abundance_data <- abundance_data[,grepl('s__',colnames(abundance_data))]
abundance_data <- abundance_data[,!grepl('t__',colnames(abundance_data))]
colnames(abundance_data) <- gsub(".*s__", "", colnames(abundance_data))

rowSums(abundance_data)
# check that they are relative abundances -> sum for sample is ~100

taxa_cut <- colnames(abundance_data)
df_for_stacked <- merge(metadata, abundance_data, by = "row.names")
row.names(df_for_stacked) <- df_for_stacked$Row.names
df_for_stacked$Row.names <- NULL

# Create a long-format dataframe where 'taxa' columns are stacked into one column
df_long <- df_for_stacked %>%
  select(all_of(taxa_cut)) %>%  # Select only the columns related to taxa
  rownames_to_column("Sample") %>%  # To keep sample names as a column
  pivot_longer(cols = -Sample ,  # Reshape from wide to long (keeping 'Sample' as identifier)
               names_to = "Taxa",  # Column name for the taxa
               values_to = "RelativeAbundance")  # Column name for the values (relative abundances)

# Calculate the mean relative abundance of each specie across all samples
mean_specie_abundance <- df_long %>%
  group_by(Taxa) %>%
  summarise(MeanAbundance = mean(RelativeAbundance)) %>%
  arrange(desc(MeanAbundance))

# Identify the 30 species with highest mean relative abundance
top_species <- mean_specie_abundance$Taxa[1:30]

# Modify df_long to keep only the 30 species with highest mean relative abundance, grouping others as "Other"
df_long <- df_long %>%
  mutate(Taxa = ifelse(Taxa %in% top_species, Taxa, "Other"))

# Set "Other" as the last factor level in Taxa
df_long$Taxa <- factor(df_long$Taxa, levels = c(top_species, "Other"))

# for each sample, compute the sum of the relative abundance of all the species classified as "Other" 
df_long <- df_long %>%
  group_by(Sample, Taxa) %>%
  summarise(RelativeAbundance = sum(RelativeAbundance)) %>%
  ungroup()

df_long$Cluster <- unlist(lapply(as.list(df_long$Sample), function(x) {metadata[x,'Cluster']}))

# order samples according to the result of hierarchical clustering
df_long$Sample <- factor(df_long$Sample, levels = c("S55101_2", "S55100_1","S55107_8","S55106_7","S55108_9","S55104_5", "S55105_6",
                                                    "S55102_3","S55103_4","S55109_10"))

n <- 31
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
colori <- col_vector
palette(colori)

# Create the bar plot with samples on the x-axis and cumulative relative abundances on the y-axis
pdf('../results_r/stacked_barplot_species.pdf', width = 8, height = 5)
ggplot(df_long, aes(x = Sample, y = RelativeAbundance, fill = Taxa)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +  # Stacked bar plot
  labs(title = "30 species with highest Mean Relative Abundance by Sample", x = "Sample", y = "Mean Relative Abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3)) +  # Rotate x-axis labels
  scale_fill_manual(values = unique(palette(colori))[1:length(unique(df_long$Taxa))]) 
dev.off()


####species welsh test####
## Differential abundance analysis on species

# Filter samples by "Cluster" to include only "C1" and "C2"
filtered_df <- df[df$Cluster %in% c("C1", "C2"), ]

# Extract the abundance data (only taxa columns) and remove taxa that are zero across all samples
abundance_data <- filtered_df[, taxa]
non_zero_taxa <- colSums(abundance_data != 0) > 0
abundance_data <- abundance_data[, non_zero_taxa]

## keep only species
abundance_data <- abundance_data[,grepl('s__',colnames(abundance_data))]
abundance_data <- abundance_data[,!grepl('t__',colnames(abundance_data))]
abundance_data <- abundance_data[,!grepl('GGB',colnames(abundance_data))]

# Update the filtered data frame to include only the non-zero taxa
filtered_df <- cbind(filtered_df[, meta], abundance_data)

# Initialize vectors to store p-values, effect sizes, and direction
p_values <- c()
effect_sizes <- c()
directions <- c()

# Perform welch t test, calculate effect size, and determine direction for each taxon
for (taxon in names(abundance_data)) {
  # Create a temporary data frame for the current taxon
  temp_df <- data.frame(
    Abundance = filtered_df[[taxon]],
    Cluster = filtered_df$Cluster
  )
  
  # Perform Wilcoxon rank-sum test
  test_result <- temp_df %>% t_test(Abundance ~ Cluster, var.equal=FALSE)
  p_values[taxon] <- test_result$p
  
  # Calculate effect size for the current taxon using wilcox_effsize
  effect_size_result <- temp_df %>% cohens_d(Abundance ~ Cluster, var.equal = FALSE)
  effect_sizes[taxon] <- effect_size_result$effsize
  
  # Determine the direction of enrichment by comparing group means
  mean_C1 <- mean(temp_df$Abundance[temp_df$Cluster == "C1"])
  mean_C2 <- mean(temp_df$Abundance[temp_df$Cluster == "C2"])
  if (mean_C1 > mean_C2) {
    directions[taxon] <- "Enriched in C1"
  } else if (mean_C2 > mean_C1) {
    directions[taxon] <- "Enriched in C2"
  } else {
    directions[taxon] <- "No Difference"
  }
}

# Adjust p-values for multiple testing to control for FDR
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Combine results into a data frame with all taxa, including effect sizes and direction
results <- data.frame(
  Taxon = names(adjusted_p_values),
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values,
  Effect_Size = effect_sizes,
  Direction = directions
)

# Extract the taxa names after "s__"
results$Taxon2 <- gsub(".*s__", "", results$Taxon)

write.table(results, file = '../results_r/welsh_test_on_species.txt', quote=F, sep='\t',row.names = F, col.names = T)

# Filter the results for p-values < 0.05
filtered_results <- results[results$P_Value <= 0.05, ]

# Extract the enriched taxa for Control and CRC groups
enriched_C1 <- filtered_results[filtered_results$Direction == "Enriched in C1", ]
enriched_C2 <- filtered_results[filtered_results$Direction == "Enriched in C2", ]

# Adjust the effect sizes for direction (positive for Control, negative for CRC)
enriched_C1$Effect_Size <- abs(enriched_C1$Effect_Size)  # Keep effect size positive for control
enriched_C2$Effect_Size <- -abs(enriched_C2$Effect_Size)  # Make effect size negative for CRC

# Prepare data for plotting
plot_data <- rbind(
  data.frame(Taxon = enriched_C1$Taxon2, 
             Effect_Size = enriched_C1$Effect_Size, 
             P_Value = enriched_C1$P_Value,
             Cluster = "C1"),
  data.frame(Taxon = enriched_C2$Taxon2, 
             Effect_Size = enriched_C2$Effect_Size, 
             P_Value = enriched_C2$P_Value,
             Cluster = "C2")
)

plot_data$P_Value <- round(plot_data$P_Value,3)

# Plot the bar plot
pdf('../results_r/effect_size_species.pdf', height = 3, width = 6)
ggplot(plot_data, aes(x = reorder(Taxon, Effect_Size), y = Effect_Size, fill = Cluster)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("C1" = "#ffb703", "C2" = "#fb8500")) +
  labs(x='',y = "Effect Size", title = "Effect Sizes of Enriched Species") +
  theme_classic() +
  geom_text(aes(label = P_Value), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(plot_data$Effect_Size >= 0, 1.2, -0.2)) +
  theme(axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))
dev.off()

####species correlation cytokines####
# correlation between all differentially abundant species (in either C1 or C2) and cytockines
cytockines <- data.frame(readxl::read_excel('../IMMUBAC-results CD45.xlsx',sheet = 'frequenze %'))
cytockines$X.CD45..IFN <- cytockines$X.CD45..IFN*100
cytockines$X.CD45..IL10 <- cytockines$X.CD45..IL10*100
cytockines$X.CD45..IL17 <- cytockines$X.CD45..IL17*100
cytockines$X.CD45..IL22 <- cytockines$X.CD45..IL22*100
cytockines$X.CD45..TNF <- cytockines$X.CD45..TNF*100

colnames(abundance_data) <- gsub(".*s__", "", colnames(abundance_data))

abundance_enriched_species_C1_C2 <- abundance_data[,filtered_results$Taxon2]
abundance_enriched_species_C1_C2_all_recipients <- do.call("rbind", replicate( 
  16, abundance_enriched_species_C1_C2, simplify = FALSE)) 

rownames(abundance_enriched_species_C1_C2_all_recipients) <- unlist(lapply(as.list(1:16), function(x) {
  paste0(rownames(abundance_enriched_species_C1_C2),'_R',as.character(x))
}))

cytockines$sample <- paste0(cytockines$donors,'_',cytockines$recipients)
rownames(cytockines) <- cytockines$sample
cytockines$sample <- NULL
cytockines$donors <- NULL
cytockines$recipients <- NULL

common <- intersect(rownames(abundance_enriched_species_C1_C2_all_recipients),rownames(cytockines))
abundance_enriched_species_C1_C2_all_recipients <- abundance_enriched_species_C1_C2_all_recipients[rownames(abundance_enriched_species_C1_C2_all_recipients) %in% common,]
cytockines <- cytockines[rownames(cytockines) %in% common,]
cytockines <- cytockines[rownames(abundance_enriched_species_C1_C2_all_recipients),]

ct1 = corr.test(abundance_enriched_species_C1_C2_all_recipients,cytockines,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p_adj<-ct1$p.adj
clinic_p_sign <- ifelse(clinic_p <= 0.05,'*','')

write.table(clinic_r, file = '../results_r/spearman_cor_species_frequency_cytockines.txt', quote=F, sep='\t',row.names = T, col.names = T)

# correlation between all differentially abundant species (in either C1 or C2) and FC cytockines
cytockines_FC <- data.frame(readxl::read_excel('../IMMUBAC-results CD45.xlsx',sheet = 'FC'))

#colnames(abundance_data) <- gsub(".*s__", "", colnames(abundance_data))

abundance_enriched_species_C1_C2 <- abundance_data[,filtered_results$Taxon2]
abundance_enriched_species_C1_C2_all_recipients <- do.call("rbind", replicate( 
  16, abundance_enriched_species_C1_C2, simplify = FALSE)) 

rownames(abundance_enriched_species_C1_C2_all_recipients) <- unlist(lapply(as.list(1:16), function(x) {
  paste0(rownames(abundance_enriched_species_C1_C2),'_R',as.character(x))
}))

cytockines_FC$sample <- paste0(cytockines_FC$donors,'_',cytockines_FC$recipients)
rownames(cytockines_FC) <- cytockines_FC$sample
cytockines_FC$sample <- NULL
cytockines_FC$donors <- NULL
cytockines_FC$recipients <- NULL

common <- intersect(rownames(abundance_enriched_species_C1_C2_all_recipients),rownames(cytockines_FC))
abundance_enriched_species_C1_C2_all_recipients <- abundance_enriched_species_C1_C2_all_recipients[rownames(abundance_enriched_species_C1_C2_all_recipients) %in% common,]
cytockines_FC <- cytockines_FC[rownames(cytockines_FC) %in% common,]
cytockines_FC <- cytockines_FC[rownames(abundance_enriched_species_C1_C2_all_recipients),]

ct1 = corr.test(abundance_enriched_species_C1_C2_all_recipients,cytockines_FC,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p_adj<-ct1$p.adj
clinic_p_sign <- ifelse(clinic_p <= 0.05,'*','')

write.table(clinic_r, file = '../results_r/spearman_cor_species_FC_cytockines.txt', quote=F, sep='\t',row.names = T, col.names = T)

# only IL22 is significantly correlated with some taxa
# Initialize vectors to store the taxa names, correlation coefficients, and p-values
taxa_names <- character(length(colnames(abundance_enriched_species_C1_C2)))
correlations <- numeric(length(colnames(abundance_enriched_species_C1_C2)))
p_values <- numeric(length(colnames(abundance_enriched_species_C1_C2)))

# Loop through each filtered taxon and compute the correlation with IL22 using Spearman correlation
for (i in 1:ncol(abundance_enriched_species_C1_C2_all_recipients)) {
  taxon <- colnames(abundance_enriched_species_C1_C2_all_recipients)[i]
  
  # Perform the Spearman correlation test between the taxon and BMI
  cor_test <- cor.test(abundance_enriched_species_C1_C2_all_recipients[,taxon], cytockines_FC$X.CD45..IL22, method = "spearman")
  
  # Store the results in the respective vectors
  taxa_names[i] <- taxon
  correlations[i] <- cor_test$estimate
  p_values[i] <- cor_test$p.value
}

# Apply FDR correction to p-values
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Create a data frame to store the results with columns for taxa, correlation, p-value, and adjusted p-value
cor_df <- data.frame(
  Taxon = taxa_names,
  Correlation = correlations,
  P_value = p_values,
  Adjusted_P_value = adjusted_p_values
)

write.table(cor_df, file = '../results_r/spearman_cor_species_FC_IL22.txt', quote=F, sep='\t',row.names = F, col.names = T)

# Order the data by absolute correlation values
cor_df$Abs_Correlation <- abs(cor_df$Correlation)

# Sort the data frame by absolute correlation values in descending order
cor_df_sorted <- cor_df[order(cor_df$Abs_Correlation, decreasing = TRUE), ]

# Select the top 20 correlations (both positive and negative)
#top_20_correlations <- cor_df_sorted[1:20, ]

# Apply a threshold of 0.2 for significance based on  p-value
cor_df_sorted_significant<- cor_df_sorted[cor_df_sorted$P_value <= 0.05,]

cor_df_sorted_significant$P_value <- round(cor_df_sorted_significant$P_value,3)

# Plot bar plot of the top 20 correlations, colored by significance
pdf('../results_r/correlation_species_FC_IL22.pdf', height = 2, width = 5)
ggplot(cor_df_sorted_significant, aes(x = reorder(Taxon, Correlation), y = Correlation)) + 
  geom_bar(stat = "identity",  fill = 'grey') +
  coord_flip() + 
  labs(title = "Correlations with FC CD45+/IL22", 
       x = "", y = "Spearman Correlation") + 
  geom_text(aes(label = P_value), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(cor_df_sorted_significant$Correlation >= 0, 1.2, -0.2)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10), title = element_text(size=10))  # Adjust text size for readability
dev.off()


####pathways beta####
pathway_table <- data.frame(read_tsv('merged_pathabundance/all_pathabundance.tsv'))
rownames(pathway_table) <- pathway_table[,1]
pathway_table[,1] <- NULL
pathway_table[,1] <- NULL
pathway_table[,1] <- NULL

colnames(pathway_table) <- rownames(metadata)[3:nrow(metadata)]
metadata <- read.table('../metadata.txt')

pathway_table <- as.data.frame(t(pathway_table))
pathways <- colnames(pathway_table)
meta <- colnames(metadata)

df_path <- merge(metadata, pathway_table, by = "row.names")

row.names(df_path) <- df_path$Row.names
df_path$Row.names <- NULL

# We only consider pathway abundance data in aggregated form, i.e., we do not consider pathways stratified by species
pathwfilt <- pathways[!grepl("\\|", pathways)] # Exclude pathways stratified by species
pathwfilt <- pathwfilt[!grepl("UNMAPPED", pathwfilt)] # Exclude "UNMAPPED"
pathwfilt <- pathwfilt[!grepl("UNINTEGRATED", pathwfilt)] # Exclude "UNINTEGRATED"

## compute beta diversity based on pathways
# Beta diversity analysis: Bray-Curtis distance calculation
beta_dist <- vegdist(select(df_path, all_of(pathwfilt)), method = "bray")

beta_dist_table <- as.matrix(beta_dist)
write.table(beta_dist_table, file = '../results_r/bray_curtis_based_on_pathways.txt', quote=F, sep='\t',row.names = T, col.names = T)

# PCoA analysis
pcoa_results <- cmdscale(beta_dist, eig = TRUE, k = 2)  # Perform PCoA with 2 components
pcoa_df <- data.frame(
  MDS1 = pcoa_results$points[, 1],
  MDS2 = pcoa_results$points[, 2],
  cluster = df_path$Cluster  # Add non_westernized for coloring
)

# Set the order of non_westernized factor
pcoa_df$cluster <- factor(pcoa_df$cluster, levels = c("C1", "C2"))

# Define color palette
colors <- c("C1" = "#ffb703", "C2" = "#fb8500")

# Perform PERMANOVA (adonis2) to compare groups based on Bray-Curtis distance
set.seed(1004)
adonis_results <- adonis2(beta_dist ~ df_path$Cluster, permutations = 999)
# 0.079

# PCoA plot with jittered points and p-value
pdf('../results_r/PCA_pathways.pdf', width = 5, height = 5)
ggplot(pcoa_df, aes(x = MDS1, y = MDS2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +  # Add points with opacity
  scale_color_manual(values = colors) +
  #stat_ellipse(aes(group = StudyCondition), level = 0.95, type = "t", size = 1) +  # Add ellipses
  labs(title = "PCoA of Beta Diversity (Bray-Curtis)", x = "MDS1", y = "MDS2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  annotate("text", x = Inf, y = Inf, label = paste("p = ", round(adonis_results$`Pr(>F)`[1], 3)),
           hjust = 2, vjust = 1.1, size = 5, color = "black")  # Add p-value from PERMANOVA
dev.off()

####pathways welsh test####
# Differential abundance analysis on pathways

# Filter samples by "Cluster" to include only "C1" and "C2"
filtered_dfPW <- df_path[df_path$Cluster %in% c("C1", "C2"), ]

# Extract the pathway abundance data (only pathwfilt columns) and remove pathways that are zero across all samples
pathw_data <- filtered_dfPW[, pathwfilt]
non_zero_pathw <- colSums(pathw_data != 0) > 0
pathw_data <- pathw_data[, non_zero_pathw]

# Update the filtered data frame to include only the non-zero taxa
filtered_dfPW <- cbind(filtered_dfPW[, meta], pathw_data)

# Initialize vectors to store p-values, effect sizes, and direction
p_values <- c()
effect_sizes <- c()
directions <- c()

# Perform welsh test, calculate effect size, and determine direction for each taxon
for (pathway in names(pathw_data)) {
  # Create a temporary data frame for the current taxon
  temp_df <- data.frame(
    Abundance = filtered_dfPW[[pathway]],
    Cluster = filtered_dfPW$Cluster
  )
  
  # Perform Wilcoxon rank-sum test
  test_result <- temp_df %>% t_test(Abundance ~ Cluster, var.equal=FALSE)
  p_values[pathway] <- test_result$p
  
  # Calculate effect size for the current taxon using wilcox_effsize
  effect_size_result <- temp_df %>% cohens_d(Abundance ~ Cluster, var.equal = FALSE)
  effect_sizes[pathway] <- effect_size_result$effsize
  
  # Determine the direction of enrichment by comparing group means
  mean_C1 <- mean(temp_df$Abundance[temp_df$Cluster == "C1"])
  mean_C2 <- mean(temp_df$Abundance[temp_df$Cluster == "C2"])
  if (mean_C1 > mean_C2) {
    directions[pathway] <- "Enriched in C1"
  } else if (mean_C2 > mean_C1) {
    directions[pathway] <- "Enriched in C2"
  } else {
    directions[pathway] <- "No Difference"
  }
}

# Adjust p-values for multiple testing to control for FDR
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Combine results into a data frame with all taxa, including effect sizes and direction
resultsPW <- data.frame(
  Pathway = names(adjusted_p_values),
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values,
  Effect_Size = effect_sizes,
  Direction = directions
)

write.table(resultsPW, file = '../results_r/welsh_test_on_pathways.txt', quote=F, sep='\t',row.names = T, col.names = T)

# Filter the results for p-values < 0.01
filtered_resultsPW <- resultsPW[resultsPW$P_Value <= 0.02, ]

# Extract the enriched taxa for C1 and C2 groups
enriched_C1 <- filtered_resultsPW[filtered_resultsPW$Direction == "Enriched in C1", ]
enriched_C2 <- filtered_resultsPW[filtered_resultsPW$Direction == "Enriched in C2", ]

# Adjust the effect sizes for direction (positive for Control, negative for CRC)
enriched_C1$Effect_Size <- abs(enriched_C1$Effect_Size)  # Keep effect size positive for no
enriched_C2$Effect_Size <- -abs(enriched_C2$Effect_Size)  # Make effect size negative for yes

# Prepare data for plotting
plot_data <- rbind(
  data.frame(Pathway = enriched_C1$Pathway,
             Effect_Size = enriched_C1$Effect_Size, 
             P_Value = enriched_C1$P_Value,
             Cluster = "C1"),
  data.frame(Pathway = enriched_C2$Pathway, 
             Effect_Size = enriched_C2$Effect_Size, 
             P_Value = enriched_C2$P_Value,
             Cluster = "C2")
)

plot_data$Pathway <- gsub('.*: ','', plot_data$Pathway)

plot_data$P_Value <- round(plot_data$P_Value,3)

# Plot the bar plot
pdf('../results_r/effect_size_pathways.pdf', height = 6, width = 8)
ggplot(plot_data, aes(x = reorder(Pathway, Effect_Size), y = Effect_Size, fill = Cluster)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = P_Value), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(plot_data$Effect_Size >= 0, 1.2, -0.2)) +
  scale_fill_manual(values = c("C1" = "#ffb703", "C2" = "#fb8500")) +
  labs(x='',y = "Effect Size", title = "Effect Sizes of Enriched Pathways") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))
dev.off()


####pathways correlation cytokines####
# correlation between all differentially abundant pathways (in either C1 or C2) and FC of cytockines
cytockines_FC <- data.frame(readxl::read_excel('../IMMUBAC-results CD45.xlsx',sheet = 'FC'))

colnames(pathw_data) <- gsub('.*: ','', colnames(pathw_data))

abundance_enriched_pathways_C1_C2 <- pathw_data[,plot_data$Pathway]
abundance_enriched_pathways_C1_C2_all_recipients <- do.call("rbind", replicate( 
  16, abundance_enriched_pathways_C1_C2, simplify = FALSE)) 

rownames(abundance_enriched_pathways_C1_C2_all_recipients) <- unlist(lapply(as.list(1:16), function(x) {
  paste0(rownames(abundance_enriched_pathways_C1_C2),'_R',as.character(x))
}))

cytockines_FC$sample <- paste0(cytockines_FC$donors,'_',cytockines_FC$recipients)
rownames(cytockines_FC) <- cytockines_FC$sample
cytockines_FC$sample <- NULL
cytockines_FC$donors <- NULL
cytockines_FC$recipients <- NULL

common <- intersect(rownames(abundance_enriched_pathways_C1_C2_all_recipients),rownames(cytockines_FC))
abundance_enriched_pathways_C1_C2_all_recipients <- abundance_enriched_pathways_C1_C2_all_recipients[rownames(abundance_enriched_pathways_C1_C2_all_recipients) %in% common,]
cytockines_FC <- cytockines_FC[rownames(cytockines_FC) %in% common,]
cytockines_FC <- cytockines_FC[rownames(abundance_enriched_pathways_C1_C2_all_recipients),]

ct1 = corr.test(abundance_enriched_pathways_C1_C2_all_recipients,cytockines_FC,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p_adj<-ct1$p.adj
clinic_p_sign <- ifelse(clinic_p <= 0.05,'*','')

write.table(clinic_r, file = '../results_r/spearman_cor_pathways_FC_cytockines.txt', quote=F, sep='\t',row.names = T, col.names = T)

# only IL22 is significantly correlated with some taxa
# Initialize vectors to store the taxa names, correlation coefficients, and p-values
path_names <- character(length(colnames(abundance_enriched_pathways_C1_C2_all_recipients)))
correlations <- numeric(length(colnames(abundance_enriched_pathways_C1_C2_all_recipients)))
p_values <- numeric(length(colnames(abundance_enriched_pathways_C1_C2_all_recipients)))

# Loop through each filtered taxon and compute the correlation with IL22 using Spearman correlation
for (i in 1:ncol(abundance_enriched_pathways_C1_C2_all_recipients)) {
  path <- colnames(abundance_enriched_pathways_C1_C2_all_recipients)[i]
  
  # Perform the Spearman correlation test between the taxon and BMI
  cor_test <- cor.test(abundance_enriched_pathways_C1_C2_all_recipients[,path], cytockines_FC$X.CD45..IL22, method = "spearman")
  
  # Store the results in the respective vectors
  path_names[i] <- path
  correlations[i] <- cor_test$estimate
  p_values[i] <- cor_test$p.value
}

# Apply FDR correction to p-values
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Create a data frame to store the results with columns for taxa, correlation, p-value, and adjusted p-value
cor_df <- data.frame(
  Pathway = path_names,
  Correlation = correlations,
  P_value = p_values,
  Adjusted_P_value = adjusted_p_values
)

write.table(cor_df, file = '../results_r/spearman_cor_pathways_FC_IL22.txt', quote=F, sep='\t',row.names = F, col.names = T)

# Order the data by absolute correlation values
cor_df$Abs_Correlation <- abs(cor_df$Correlation)

# Sort the data frame by absolute correlation values in descending order
cor_df_sorted <- cor_df[order(cor_df$Abs_Correlation, decreasing = TRUE), ]

# Select the top 20 correlations (both positive and negative)
#top_20_correlations <- cor_df_sorted[1:20, ]

# Apply a threshold of 0.2 for significance based on  p-value
cor_df_sorted_significant <- cor_df_sorted[cor_df_sorted$Adjusted_P_value <= 0.05,]

cor_df_sorted_significant$Adjusted_P_value <- round(cor_df_sorted_significant$Adjusted_P_value,3)

# Plot bar plot of the top 20 correlations, colored by significance
pdf('../results_r/correlation_pathways_FC_IL22.pdf', height = 4, width = 7.5)
ggplot(cor_df_sorted_significant, aes(x = reorder(Pathway, Correlation), y = Correlation)) + 
  geom_bar(stat = "identity",  fill='gray') +
  coord_flip() + 
  labs(title = "Correlations with FC CD45+/IL22", 
       x = "", y = "Spearman Correlation") + 
  geom_text(aes(label = Adjusted_P_value), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(cor_df_sorted_significant$Correlation >= 0, 1.2, -0.2)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10), title = element_text(size=10))  # Adjust text size for readability
dev.off()

####gene families beta####
gene_fam_table <- data.frame(read_tsv('merged_genefamilies_KO_renamed/all_genefamilies_KO_renamed.tsv'))
rownames(gene_fam_table) <- gene_fam_table[,1]
gene_fam_table[,1] <- NULL
gene_fam_table[,1] <- NULL
gene_fam_table[,1] <- NULL

colnames(gene_fam_table) <- rownames(metadata)[3:nrow(metadata)]

gene_fam_table <- as.data.frame(t(gene_fam_table))
gene_families <- colnames(gene_fam_table)
meta <- colnames(metadata)

df_gene_fam <- merge(metadata, gene_fam_table, by = "row.names")

row.names(df_gene_fam) <- df_gene_fam$Row.names
df_gene_fam$Row.names <- NULL

gene_fami_filt <- gene_families[!grepl("\\|", gene_families)] # Exclude gene families stratified by species
gene_fami_filt <- gene_fami_filt[!grepl("UNMAPPED", gene_fami_filt)] # Exclude "UNMAPPED"
gene_fami_filt <- gene_fami_filt[!grepl("UNINTEGRATED", gene_fami_filt)] # Exclude "UNINTEGRATED"

## compute beta diversity based on gene families
# Beta diversity analysis: Bray-Curtis distance calculation
beta_dist <- vegdist(select(df_gene_fam, all_of(gene_fami_filt)), method = "bray")

beta_dist_table <- as.matrix(beta_dist)
write.table(beta_dist_table, file = '../results_r/bray_curtis_based_on_gene_families.txt', quote=F, sep='\t',row.names = T, col.names = T)

# PCoA analysis
pcoa_results <- cmdscale(beta_dist, eig = TRUE, k = 2)  # Perform PCoA with 2 components
pcoa_df <- data.frame(
  MDS1 = pcoa_results$points[, 1],
  MDS2 = pcoa_results$points[, 2],
  cluster = df_gene_fam$Cluster  # Add non_westernized for coloring
)

# Set the order of non_westernized factor
pcoa_df$cluster <- factor(pcoa_df$cluster, levels = c("C1", "C2"))

# Define color palette
colors <- c("C1" = "#ffb703", "C2" = "#fb8500")

# Perform PERMANOVA (adonis2) to compare groups based on Bray-Curtis distance
set.seed(1004)
adonis_results <- adonis2(beta_dist ~ df_gene_fam$Cluster, permutations = 999)
# 0.182

# PCoA plot with jittered points and p-value
pdf('../results_r/PCA_gene_families.pdf', width = 5, height = 5)
ggplot(pcoa_df, aes(x = MDS1, y = MDS2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +  # Add points with opacity
  scale_color_manual(values = colors) +
  #stat_ellipse(aes(group = StudyCondition), level = 0.95, type = "t", size = 1) +  # Add ellipses
  labs(title = "PCoA of Beta Diversity (Bray-Curtis)", x = "MDS1", y = "MDS2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  annotate("text", x = Inf, y = Inf, label = paste("p = ", round(adonis_results$`Pr(>F)`[1], 3)),
           hjust = 3.6, vjust = 1.1, size = 5, color = "black")  # Add p-value from PERMANOVA
dev.off()

####gene families welsh test####
# Differential abundance analysis on pathways: the Wilcoxon rank-sum test

# Filter samples by "Cluster" to include only "C1" and "C2"
filtered_gene_fam <- df_gene_fam[df_gene_fam$Cluster %in% c("C1", "C2"), ]

# Extract the gene fam. abundance data (only gene_fami_filt columns) and remove gene fam. that are zero across all samples
gene_fam_data <- filtered_gene_fam[, gene_fami_filt]
non_zero_gene_fam <- colSums(gene_fam_data != 0) > 0
gene_fam_data <- gene_fam_data[, non_zero_gene_fam]

# Update the filtered data frame to include only the non-zero taxa
filtered_gene_fam <- cbind(filtered_gene_fam[, meta], gene_fam_data)

# Initialize vectors to store p-values, effect sizes, and direction
p_values <- c()
effect_sizes <- c()
directions <- c()

# Perform Wilcoxon rank-sum test, calculate effect size, and determine direction for each taxon
for (fam in names(gene_fam_data)) {
  # Create a temporary data frame for the current taxon
  temp_df <- data.frame(
    Abundance = filtered_gene_fam[[fam]],
    Cluster = filtered_gene_fam$Cluster
  )
  
  # Perform Wilcoxon rank-sum test
  test_result <- temp_df %>% t_test(Abundance ~ Cluster, var.equal=FALSE)
  p_values[fam] <- test_result$p
  
  # Calculate effect size for the current taxon using wilcox_effsize
  effect_size_result <- temp_df %>% cohens_d(Abundance ~ Cluster, var.equal = FALSE)
  effect_sizes[fam] <- effect_size_result$effsize
  
  # Determine the direction of enrichment by comparing group means
  mean_C1 <- mean(temp_df$Abundance[temp_df$Cluster == "C1"])
  mean_C2 <- mean(temp_df$Abundance[temp_df$Cluster == "C2"])
  if (mean_C1 > mean_C2) {
    directions[fam] <- "Enriched in C1"
  } else if (mean_C2 > mean_C1) {
    directions[fam] <- "Enriched in C2"
  } else {
    directions[fam] <- "No Difference"
  }
}

# Adjust p-values for multiple testing to control for FDR
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Combine results into a data frame with all taxa, including effect sizes and direction
results_gene_fam <- data.frame(
  Gene_family = names(adjusted_p_values),
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values,
  Effect_Size = effect_sizes,
  Direction = directions
)

write.table(results_gene_fam, file = '../results_r/welsh_test_on_gene_families.txt', quote=F, sep='\t',row.names = T, col.names = T)

# Filter the results for p-values < 0.02
filtered_results_gene_fam <- results_gene_fam[results_gene_fam$P_Value <= 0.01, ]

# Extract the enriched taxa for C1 and C2 groups
enriched_C1 <- filtered_results_gene_fam[filtered_results_gene_fam$Direction == "Enriched in C1", ]
enriched_C2 <- filtered_results_gene_fam[filtered_results_gene_fam$Direction == "Enriched in C2", ]

# Adjust the effect sizes for direction (positive for Control, negative for CRC)
enriched_C1$Effect_Size <- abs(enriched_C1$Effect_Size)  # Keep effect size positive for no
enriched_C2$Effect_Size <- -abs(enriched_C2$Effect_Size)  # Make effect size negative for yes

# Prepare data for plotting
plot_data <- rbind(
  data.frame(Gene_family = enriched_C1$Gene_family,
             Effect_Size = enriched_C1$Effect_Size, 
             P_Value = enriched_C1$P_Value,
             Cluster = "C1"),
   data.frame(Gene_family = enriched_C2$Gene_family, 
              Effect_Size = enriched_C2$Effect_Size, 
              P_Value = enriched_C2$P_Value,
             Cluster = "C2")
)

plot_data$Gene_family <- gsub('\\[.*','', plot_data$Gene_family)
plot_data$Gene_family <- gsub(': NO_NAME','', plot_data$Gene_family)

plot_data$P_Value <- round(plot_data$P_Value,3)

# Plot the bar plot
pdf('../results_r/effect_size_gene_families.pdf', height = 10, width = 10)
ggplot(plot_data, aes(x = reorder(Gene_family, Effect_Size), y = Effect_Size, fill = Cluster)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = P_Value), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(plot_data$Effect_Size >= 0, 1.2, -0.2)) +
  scale_fill_manual(values = c("C1" = "#ffb703", "C2" = "#fb8500")) +
  labs(x='',y = "Effect Size", title = "Effect Sizes of Enriched Gene families") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))
dev.off()


####gene families correlation cytokines####
# correlation between all differentially abundant pathways (in either C1 or C2) and cytockines
cytockines_FC <- data.frame(readxl::read_excel('../IMMUBAC-results CD45.xlsx',sheet = 'FC'))

colnames(gene_fam_data) <- gsub('\\[.*','', colnames(gene_fam_data))
colnames(gene_fam_data)  <- gsub(': NO_NAME','', colnames(gene_fam_data))

abundance_enriched_gene_fam_C1_C2 <- gene_fam_data[,plot_data$Gene_family]
abundance_enriched_gene_fam_C1_C2_all_recipients <- do.call("rbind", replicate( 
  16, abundance_enriched_gene_fam_C1_C2, simplify = FALSE)) 

rownames(abundance_enriched_gene_fam_C1_C2_all_recipients) <- unlist(lapply(as.list(1:16), function(x) {
  paste0(rownames(abundance_enriched_gene_fam_C1_C2),'_R',as.character(x))
}))

cytockines_FC$sample <- paste0(cytockines_FC$donors,'_',cytockines_FC$recipients)
rownames(cytockines_FC) <- cytockines_FC$sample
cytockines_FC$sample <- NULL
cytockines_FC$donors <- NULL
cytockines_FC$recipients <- NULL

common <- intersect(rownames(abundance_enriched_gene_fam_C1_C2_all_recipients),rownames(cytockines_FC))
abundance_enriched_gene_fam_C1_C2_all_recipients <- abundance_enriched_gene_fam_C1_C2_all_recipients[rownames(abundance_enriched_gene_fam_C1_C2_all_recipients) %in% common,]
cytockines_FC <- cytockines_FC[rownames(cytockines_FC) %in% common,]
cytockines_FC <- cytockines_FC[rownames(abundance_enriched_gene_fam_C1_C2_all_recipients),]

ct1 = corr.test(abundance_enriched_gene_fam_C1_C2_all_recipients,cytockines_FC,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p_adj<-ct1$p.adj
clinic_p_sign <- ifelse(clinic_p <= 0.07,'*','')

write.table(clinic_r, file = '../results_r/spearman_cor_gene_families_FC_cytockines.txt', quote=F, sep='\t',row.names = T, col.names = T)

# only IL22 is significantly correlated with some taxa
# Initialize vectors to store the taxa names, correlation coefficients, and p-values
gene_fam_names <- character(length(colnames(abundance_enriched_gene_fam_C1_C2_all_recipients)))
correlations <- numeric(length(colnames(abundance_enriched_gene_fam_C1_C2_all_recipients)))
p_values <- numeric(length(colnames(abundance_enriched_gene_fam_C1_C2_all_recipients)))

# Loop through each filtered taxon and compute the correlation with IL22 using Spearman correlation
for (i in 1:ncol(abundance_enriched_gene_fam_C1_C2_all_recipients)) {
  gene_fam <- colnames(abundance_enriched_gene_fam_C1_C2_all_recipients)[i]
  
  # Perform the Spearman correlation test between the taxon and BMI
  cor_test <- cor.test(abundance_enriched_gene_fam_C1_C2_all_recipients[,gene_fam], cytockines_FC$X.CD45..IL22, method = "spearman")
  
  # Store the results in the respective vectors
  gene_fam_names[i] <- gene_fam
  correlations[i] <- cor_test$estimate
  p_values[i] <- cor_test$p.value
}

# Apply FDR correction to p-values
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Create a data frame to store the results with columns for taxa, correlation, p-value, and adjusted p-value
cor_df <- data.frame(
  Gene_families = gene_fam_names,
  Correlation = correlations,
  P_value = p_values,
  Adjusted_P_value = adjusted_p_values
)

write.table(cor_df, file = '../results_r/spearman_cor_gene_families_FC_IL22.txt', quote=F, sep='\t',row.names = F, col.names = T)

# Order the data by absolute correlation values
cor_df$Abs_Correlation <- abs(cor_df$Correlation)

# Sort the data frame by absolute correlation values in descending order
cor_df_sorted <- cor_df[order(cor_df$Abs_Correlation, decreasing = TRUE), ]

# Select the top 20 correlations (both positive and negative)
#top_20_correlations <- cor_df_sorted[1:20, ]

# Apply a threshold of 0.2 for significance based on  p-value
cor_df_sorted_significant <- cor_df_sorted[cor_df_sorted$Adjusted_P_value <= 0.05, ]

cor_df_sorted_significant$Adjusted_P_value <- round(cor_df_sorted_significant$Adjusted_P_value,3)

# Plot bar plot of the top 20 correlations, colored by significance
pdf('../results_r/correlation_gene_fam_FC_IL22.pdf', height = 10, width = 8)
ggplot(cor_df_sorted_significant, aes(x = reorder(Gene_families, Correlation), y = Correlation)) + 
  geom_bar(stat = "identity", color = "black", fill='gray') +
  coord_flip() + 
  labs(title = "Correlations with FC CD45+/IL22", 
       x = "", y = "Spearman Correlation") + 
  geom_text(aes(label = Adjusted_P_value), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(cor_df_sorted_significant$Correlation >= 0, 1.2, -0.2)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), title = element_text(size=10))  # Adjust text size for readability
dev.off()


####GMM welsh test####
## GMM modules enriched between C1 and C2
gene_fam_table <- data.frame(read_tsv('merged_genefamilies_KO/all_genefamilies_KO.tsv'))
rownames(gene_fam_table) <- gene_fam_table[,1]
gene_fam_table[,1] <- NULL
gene_fam_table[,1] <- NULL
gene_fam_table[,1] <- NULL

colnames(gene_fam_table) <- rownames(metadata)[3:nrow(metadata)]

gene_families <- rownames(gene_fam_table)

gene_fami_filt <- gene_families[!grepl("\\|", gene_families)] # Exclude gene families stratified by species
gene_fami_filt <- gene_fami_filt[!grepl("UNMAPPED", gene_fami_filt)] # Exclude "UNMAPPED"
gene_fami_filt <- gene_fami_filt[!grepl("UNGROUPED", gene_fami_filt)] # Exclude "UNINTEGRATED"

gene_fam_table <- gene_fam_table[gene_fami_filt,]

write.table(gene_fam_table, file='../metaphlan_humann/merged_genefamilies_KO/all_genefamilies_KO_for_omixer.tsv', quote=F, col.names =TRUE, row.names = TRUE, sep = '\t')

df_modules <- data.frame(Module = NA, Value=NA, Sample = NA)
for (f in list.files('./module-prediction-1741180414663/', full.names = TRUE)) {
  df_sample <- cbind(read.table(f, header = TRUE, sep = '\t')[,c(1,2)], Sample= rep(strsplit(f, split = '/')[[1]][4],nrow(read.table(f, header = TRUE, sep = '\t')[,c(1,2)])))
  df_modules <- rbind(df_modules, df_sample)
}

df_modules <- df_modules[-1,]

df_modules$Sample <- gsub('.modules','',df_modules$Sample)
df_modules$cluster <- unlist(lapply(as.list(df_modules$Sample), function(x) metadata[x,2]))

all_modules <- unique(df_modules$Module)
all_samples <- unique(df_modules$Sample)

m <- matrix(nrow=length(all_samples), ncol =length(all_modules),0)
rownames(m) <- all_samples
colnames(m) <- all_modules

for (s in rownames(m)) {
  for (mod in colnames(m)) {
    if (length(df_modules[(df_modules$Module == mod & df_modules$Sample == s),'Value']) != 0) {
    m[s,mod] <- df_modules[df_modules$Module ==mod & df_modules$Sample == s, 'Value']
    }
  }
}

modules <- colnames(m)
meta <- colnames(metadata)

df_modules <- merge(metadata, m, by = "row.names")

row.names(df_modules) <- df_modules$Row.names
df_modules$Row.names <- NULL


# Differential abundance analysis on modules: the welsh test

# Filter samples by "Cluster" to include only "C1" and "C2"
filtered_modules <- df_modules[df_modules$Cluster %in% c("C1", "C2"), ]

# Extract the modules abundance data (only gene_fami_filt columns) and remove modules that are zero across all samples
modules_data <- filtered_modules[, modules]
non_zero_modules <- colSums(modules_data != 0) > 0
modules_data <- modules_data[, non_zero_modules]

GMM <- readxl::read_excel('./../GMM_module_kegg.xlsx')

colnames(modules_data) <- unlist(lapply(as.list(as.vector(colnames(modules_data))), function(x) 
{paste(x, '-', unique(GMM[GMM$module ==x, 3]))}))

# Update the filtered data frame to include only the non-zero taxa
filtered_modules <- cbind(filtered_modules[, meta], modules_data)

# Initialize vectors to store p-values, effect sizes, and direction
p_values <- c()
effect_sizes <- c()
directions <- c()

# Perform Wilcoxon rank-sum test, calculate effect size, and determine direction for each taxon
for (mod in names(modules_data)) {
  # Create a temporary data frame for the current taxon
  temp_df <- data.frame(
    Abundance = filtered_modules[[mod]],
    Cluster = filtered_modules$Cluster
  )
  
  # Perform Wilcoxon rank-sum test
  test_result <- temp_df %>% t_test(Abundance ~ Cluster, var.equal=FALSE)
  p_values[mod] <- test_result$p
  
  # Calculate effect size for the current taxon using wilcox_effsize
  effect_size_result <- temp_df %>% cohens_d(Abundance ~ Cluster, var.equal = FALSE)
  effect_sizes[mod] <- effect_size_result$effsize
  
  # Determine the direction of enrichment by comparing group means
  mean_C1 <- mean(temp_df$Abundance[temp_df$Cluster == "C1"])
  mean_C2 <- mean(temp_df$Abundance[temp_df$Cluster == "C2"])
  if (mean_C1 > mean_C2) {
    directions[mod] <- "Enriched in C1"
  } else if (mean_C2 > mean_C1) {
    directions[mod] <- "Enriched in C2"
  } else {
    directions[mod] <- "No Difference"
  }
}

# Adjust p-values for multiple testing to control for FDR
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Combine results into a data frame with all taxa, including effect sizes and direction
results_modules <- data.frame(
  Modules = names(adjusted_p_values),
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values,
  Effect_Size = effect_sizes,
  Direction = directions
)

write.table(results_modules, file = '../results_r/welsh_test_on_GMM.txt', quote=F, sep='\t',row.names = F, col.names = T)

# Filter the results for p-values < 0.02
filtered_results_modules <- results_modules[results_modules$P_Value <= 0.05, ]

# Extract the enriched taxa for C1 and C2 groups
enriched_C1 <- filtered_results_modules[filtered_results_modules$Direction == "Enriched in C1", ]
enriched_C2 <- filtered_results_modules[filtered_results_modules$Direction == "Enriched in C2", ]

# Adjust the effect sizes for direction (positive for Control, negative for CRC)
enriched_C1$Effect_Size <- abs(enriched_C1$Effect_Size)  # Keep effect size positive for no
enriched_C2$Effect_Size <- -abs(enriched_C2$Effect_Size)  # Make effect size negative for yes

# Prepare data for plotting
plot_data <- rbind(
  data.frame(Modules = enriched_C1$Modules,
             Effect_Size = enriched_C1$Effect_Size, 
             P_Value = enriched_C1$P_Value,
             Cluster = "C1"),
  data.frame(Modules = enriched_C2$Modules,
             Effect_Size = enriched_C2$Effect_Size,
             P_Value = enriched_C2$P_Value,
             Cluster = "C2")
)

plot_data$P_Value <- round(plot_data$P_Value,3)

# Plot the bar plot
pdf('../results_r/effect_size_modules.pdf', height = 4, width = 7.5)
ggplot(plot_data, aes(x = reorder(Modules, Effect_Size), y = Effect_Size, fill = Cluster)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = P_Value), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(plot_data$Effect_Size >= 0, 1.2, -0.2)) +
  scale_fill_manual(values = c("C1" = "#ffb703", "C2" = "#fb8500")) +
  labs(x='',y = "Effect Size", title = "Effect Sizes of Enriched GMM") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))
dev.off()


####GMM correlation cytokines####
# correlation between all differentially abundant GMM (in either C1 or C2) and FC cytockines
cytockines_FC <- data.frame(readxl::read_excel('../IMMUBAC-results CD45.xlsx',sheet = 'FC'))

abundance_enriched_modules_C1_C2 <- modules_data[,plot_data$Modules]
abundance_enriched_modules_C1_C2_all_recipients <- do.call("rbind", replicate( 
  16, abundance_enriched_modules_C1_C2, simplify = FALSE)) 

rownames(abundance_enriched_modules_C1_C2_all_recipients) <- unlist(lapply(as.list(1:16), function(x) {
  paste0(rownames(abundance_enriched_modules_C1_C2),'_R',as.character(x))
}))

cytockines_FC$sample <- paste0(cytockines_FC$donors,'_',cytockines_FC$recipients)
rownames(cytockines_FC) <- cytockines_FC$sample
cytockines_FC$sample <- NULL
cytockines_FC$donors <- NULL
cytockines_FC$recipients <- NULL

common <- intersect(rownames(abundance_enriched_modules_C1_C2_all_recipients),rownames(cytockines_FC))
abundance_enriched_modules_C1_C2_all_recipients <- abundance_enriched_modules_C1_C2_all_recipients[rownames(abundance_enriched_modules_C1_C2_all_recipients) %in% common,]
cytockines_FC <- cytockines_FC[rownames(cytockines_FC) %in% common,]
cytockines_FC <- cytockines_FC[rownames(abundance_enriched_modules_C1_C2_all_recipients),]

ct1 = corr.test(abundance_enriched_modules_C1_C2_all_recipients,cytockines_FC,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p_adj<-ct1$p.adj
clinic_p_sign <- ifelse(clinic_p <= 0.05,'*','')

write.table(clinic_r, file = '../results_r/spearman_cor_modules_FC_cytockines.txt', quote=F, sep='\t',row.names = T, col.names = T)

# only IL22 is significantly correlated with some taxa
# Initialize vectors to store the taxa names, correlation coefficients, and p-values
modules_names <- character(length(colnames(abundance_enriched_modules_C1_C2_all_recipients)))
correlations <- numeric(length(colnames(abundance_enriched_modules_C1_C2_all_recipients)))
p_values <- numeric(length(colnames(abundance_enriched_modules_C1_C2_all_recipients)))

# Loop through each filtered taxon and compute the correlation with IL22 using Spearman correlation
for (i in 1:ncol(abundance_enriched_modules_C1_C2_all_recipients)) {
  mod <- colnames(abundance_enriched_modules_C1_C2_all_recipients)[i]
  
  # Perform the Spearman correlation test between the taxon and BMI
  cor_test <- cor.test(abundance_enriched_modules_C1_C2_all_recipients[,mod], cytockines_FC$X.CD45..IL22, method = "spearman")
  
  # Store the results in the respective vectors
  modules_names[i] <- mod
  correlations[i] <- cor_test$estimate
  p_values[i] <- cor_test$p.value
}

# Apply FDR correction to p-values
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Create a data frame to store the results with columns for taxa, correlation, p-value, and adjusted p-value
cor_df <- data.frame(
  Modules = modules_names,
  Correlation = correlations,
  P_value = p_values,
  Adjusted_P_value = adjusted_p_values
)

write.table(cor_df, file = '../results_r/spearman_cor_modules_FC_IL22.txt', quote=F, sep='\t',row.names = F, col.names = T)

# Order the data by absolute correlation values
cor_df$Abs_Correlation <- abs(cor_df$Correlation)

# Sort the data frame by absolute correlation values in descending order
cor_df_sorted <- cor_df[order(cor_df$Abs_Correlation, decreasing = TRUE), ]

# Select the top 20 correlations (both positive and negative)
#top_20_correlations <- cor_df_sorted[1:20, ]

# Apply a threshold of 0.2 for significance based on  p-value
cor_df_sorted_significant <- cor_df_sorted[cor_df_sorted$Adjusted_P_value <= 0.05,]

cor_df_sorted_significant$Adjusted_P_value <- round(cor_df_sorted_significant$Adjusted_P_value,3)

# Plot bar plot of the top 20 correlations, colored by significance
pdf('../results_r/correlation_modules_FC_IL22.pdf', height = 4, width = 7)
ggplot(cor_df_sorted_significant, aes(x = reorder(Modules, Correlation), y = Correlation)) + 
  geom_bar(stat = "identity", color = "black", fill='gray') +
  coord_flip() + 
  labs(title = "Correlations with FC CD45+/IL22", 
       x = "", y = "Spearman Correlation") + 
  geom_text(aes(label = Adjusted_P_value), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(cor_df_sorted_significant$Correlation >= 0, 1.2, -0.2)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), title = element_text(size=10))  # Adjust text size for readability
dev.off()


