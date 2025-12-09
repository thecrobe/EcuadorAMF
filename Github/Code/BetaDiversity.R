library(igraph)
library(microViz)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)


##### CCA Analysis #####
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)

# Prepare data
# Convert the OTU table to a matrix and clean data
load("Data/Phyloseq/AMF_Ecuador_SSU.Rdata")
amf1<-speedyseq::tax_glom(amf,taxrank = "Genus",NArm = FALSE)
otu_mat <- as.matrix(otu_table(amf))
otu_mat <- otu_mat[rowSums(otu_mat) > 0, ]  # Remove zero-sum rows
otu_mat[is.na(otu_mat)] <- 0  # Replace NA values with zero

# Extract sample data
sample_data <- data.frame(sample_data(amf))

# Read in additional environmental data
layers <- read.csv("Data/EcuadorLayers_GEE.csv", header = TRUE)
sample_data <- left_join(sample_data, layers, by = "Site_ID")  # Make sure 'SiteID' is the correct key

# Rename and adjust environmental variables
sample_data <- sample_data %>% 
  rename(MAP = CHELSA_BIO_Annual_Precipitation,
         MAT = CHELSA_BIO_Annual_Mean_Temperature) %>%
  mutate(MAT = MAT * 0.1,
         MAP = MAP * 0.1,
         P_percent = (P_mgL / 1000) * 10,
         CN = Carbon_totalOrganic_percent / (Nitrogen_percent * 10000),
         NP = (Nitrogen_percent * 10000) / P_percent)

# Statistical summaries and PCA
selected_data <- sample_data %>% select(CN, NP, pH)
pca_result <- prcomp(selected_data, center = TRUE, scale. = TRUE)
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Encoding data and preparing for models
sample_data <- sample_data %>%
  mutate(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

# Run CCA
sample_data$Landuse<-sample_data$FarmNative
cca_results <- cca(t(otu_mat) ~ MAT + MAP + Landuse + Crop + CN + NP + pH, data = sample_data)
eigenvalues<-cca_results$CCA$eig
variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Extract scores
site_scores <- as.data.frame(scores(cca_results, display = "sites"))
biplot_scores <- as.data.frame(scores(cca_results, display = "bp"))

# Merge FarmNative data for coloring
site_scores$Landuse<-sample_data$FarmNative
site_scores$Crop<-sample_data$Crop

site_scores <-  site_scores %>%
  mutate(
    Landuse_Crop = paste(sample_data$FarmNative, sample_data$Crop, sep = "_")
  )

site_scores<-site_scores %>%
  mutate(Landuse = recode(Landuse, "Farm" = "Monoculture", "Native" = "Uncultivated"))

# Create CCA plot with ggplot2 focusing on centroids, ellipses, and arrows
cca_plot<-ggplot() +
  geom_point(data = site_scores, aes(x = CCA1, y = CCA2, color = Landuse, shape = Crop), size = 2) +
  stat_summary(data = site_scores, fun = mean, geom = "point", size = 2, aes(x = CCA1, y = CCA2, color = Landuse, shape = Crop), stroke = 1.5) +
  scale_color_manual(values = c("Monoculture" = "#E93925", "Uncultivated" = "#2CB7D2")) +
  scale_shape_manual(values = c("Potato" = 16, "Maize" = 17)) +  # Customize shapes for each crop
  stat_ellipse(data = site_scores, aes(x = CCA1, y = CCA2, fill = Landuse, group = Landuse), level = 0.95, geom = "polygon", alpha = 0.2) +
  geom_segment(data = biplot_scores, aes(x = 0, y = 0, xend = CCA1*2, yend = CCA2*2), 
               arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "black", size = 0.5) +
  labs(title = "", x = "CCA Axis 1 (23.3%)", y = "CCA Axis 2 (18.8%)") +
#   geom_text(data = biplot_scores, aes(x = CCA1*2, y = CCA2*2, label = rownames(biplot_scores)), hjust = 1, vjust = 1, color = "black") +
  theme_minimal()

ggsave(plot = cca_plot,"CCA_plot.pdf",scale = 0.5)

# Assuming 'cca_results' is already created from your CCA analysis
# Permutation test for the overall model
overall_test <- anova(cca_results, permutations = 999)  # Using 999 permutations
print(overall_test)

# Permutation test for individual variables
individual_tests <- anova(cca_results, by = "terms", permutations = 999)
print(individual_tests)
varexp_covs<-0.4063+0.4520+0.3722+0.3960+0.3130+0.4230+0.4076
residual<-13.9
variance<-varexp_covs / (varexp_covs+residual)

###### Summary statistics #####
merged<-ps_filtered %>% phyloseq::merge_samples(group = "FarmNative_Crop")
merged_ra <- transform_sample_counts(merged, function(x) x / sum(x))
merged_ra_melt<-psmelt(ps_filtered)  
merged_ra_melt %>% 
  group_by(FarmNative_Crop, Species) %>%
  summarize(mean_abundance = sum(Abundance), .groups = 'drop')

genus_totals <- merged_ra_melt %>%
  group_by(FarmNative_Crop, Species) %>%
  summarize(total_abundance = sum(Abundance), .groups = 'drop')

farm_totals <- genus_totals %>%
  group_by(FarmNative_Crop) %>%
  summarize(group_total_abundance = sum(total_abundance), .groups = 'drop')

relative_abundance <- genus_totals %>%
  left_join(farm_totals, by = "FarmNative_Crop") %>%
  mutate(relative_abundance = (total_abundance / group_total_abundance) * 100) %>%
  select(FarmNative_Crop, Species, relative_abundance)

print(relative_abundance)

top_5_abundance <- relative_abundance %>%
  arrange(desc(relative_abundance)) %>%
  group_by(FarmNative_Crop) %>%
  slice_max(order_by = relative_abundance, n = 5) %>%
  ungroup()

