# -------------------------------------------------
# Step 1: Load Libraries and Data
# -------------------------------------------------
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)

# Load the AMF Dataset
load(file = "/Data/Phyloseq/AMF_Ecuador_SSU.Rdata")

# Add combined factor for FarmNative and Crop
amf@sam_data$FarmNative_Crop <- paste(amf@sam_data$FarmNative, amf@sam_data$Crop)
amf@sam_data$FarmNative_Crop <- factor(amf@sam_data$FarmNative_Crop, levels = c("Farm Maize", "Native Maize", "Farm Potato", "Native Potato"))

# Collapse taxa to genus level
amf_genus <- speedyseq::tax_glom(amf, taxrank = "Genus", NArm = FALSE)

# Extract OTU table as a matrix
otu_mat <- as.matrix(otu_table(amf_genus))
otu_mat <- otu_mat[rowSums(otu_mat) > 0, ]  # Remove zero-sum rows
otu_mat[is.na(otu_mat)] <- 0  # Replace NA values with zeros

# Extract sample metadata
sample_data <- data.frame(sample_data(amf_genus))

# Read in additional environmental layers
layers <- read.csv("Data/EcuadorLayers_GEE.csv", header = TRUE)

# Join sample metadata with environmental layers
sample_data <- left_join(sample_data, layers, by = "Site_ID")

# -------------------------------------------------
# Step 2: Perform PCA on Soil Variables
# -------------------------------------------------
# Adjust environmental variables
sample_data <- sample_data %>% 
  rename(
    MAP = CHELSA_BIO_Annual_Precipitation,
    MAT = CHELSA_BIO_Annual_Mean_Temperature
  ) %>%
  mutate(
    MAT = MAT * 0.1,  # Convert temperature to degrees Celsius
    MAP = MAP * 0.1,  # Adjust precipitation scale
    P_percent = (P_mgL / 1000) * 10,  # Convert phosphorus concentration
    CN = Carbon_totalOrganic_percent / (Nitrogen_percent * 10000),  # Carbon-to-nitrogen ratio
    NP = (Nitrogen_percent * 10000) / P_percent  # Nitrogen-to-phosphorus ratio
  )

# Select soil variables for PCA
soil_vars <- sample_data %>% dplyr::select(CN, NP, pH, Carbon_totalOrganic_percent)

# Perform PCA
soil_pca <- prcomp(soil_vars, center = TRUE, scale. = TRUE)

# Add PC1 to sample_data
sample_data <- sample_data %>%
  mutate(PC1 = soil_pca$x[, 1])

# -------------------------------------------------
# Step 3: Calculate Bray-Curtis Dissimilarity Matrix
# -------------------------------------------------
# Compute Bray-Curtis dissimilarity matrix
bray_dist <- vegdist(t(otu_mat), method = "bray")

# -------------------------------------------------
# Step 4: CAP Analysis
# -------------------------------------------------
# Run CAP analysis using PC1 and other environmental variables
capscale_results <- capscale(bray_dist ~ PC1 + MAT + MAP + FarmNative + Crop, data = sample_data)

# Summarize CAP results
cat("CAP Analysis Summary:\n")
summary(capscale_results)

# -------------------------------------------------
# Step 5: Perform PERMANOVA
# -------------------------------------------------
# Perform PERMANOVA on Bray-Curtis dissimilarity matrix
permanova_results <- adonis2(bray_dist ~ PC1 + MAT + FarmNative_Crop + Site_ID, stata=sample_data$Site_ID, 
                             data = sample_data, 
                             permutations = 10000, 
                             by = "terms")

# Print PERMANOVA results
cat("PERMANOVA Results for Each Term:\n")
print(permanova_results)

# -------------------------------------------------
# Step 6: Perform Envfit Analysis
# -------------------------------------------------
# Perform envfit analysis on CAP results
envfit_results <- envfit(capscale_results, 
                         sample_data[, c("PC1", "MAT", "FarmNative", "Crop","Site_ID")], 
                         permutations = 999)

# Extract environmental vector scores for plotting
biplot_scores <- as.data.frame(scores(envfit_results, display = "vectors")) %>%
  tibble::rownames_to_column(var = "Variable") 
# -------------------------------------------------
# Step 7: Visualization with Centroids, SDs, and Envfit Vectors
# -------------------------------------------------
# Extract site scores from CAP results
site_scores <- as.data.frame(scores(capscale_results, display = "sites"))
site_scores$FarmNative_Crop <- sample_data$FarmNative_Crop

# Calculate centroids for each FarmNative_Crop group
centroids <- site_scores %>%
  group_by(FarmNative_Crop) %>%
  summarize(CAP1_mean = mean(CAP1), CAP2_mean = mean(CAP2), .groups = "drop")

# Calculate standard deviations for each group
sd_regions <- site_scores %>%
  group_by(FarmNative_Crop) %>%
  summarize(
    CAP1_mean = mean(CAP1),
    CAP2_mean = mean(CAP2),
    CAP1_sd = sd(CAP1),
    CAP2_sd = sd(CAP2),
    .groups = "drop"
  )

# Define custom colors for groups
group_colors <- c(
  "Farm Maize" = "#E93925", "Native Maize" = "#2CB7D2",
  "Farm Potato" = "#E97D25", "Native Potato" = "#2C8BD2"
)

# Create CAP ordination plot
cap_plot<-ggplot() +
  # Background points for all sites
  geom_point(data = site_scores, aes(x = CAP1, y = CAP2, color = FarmNative_Crop), size = 2, alpha = 0.3) +
  # Add SD regions for each group
  geom_errorbar(data = sd_regions, aes(x = CAP1_mean, ymin = CAP2_mean - CAP2_sd, ymax = CAP2_mean + CAP2_sd, color = FarmNative_Crop), width = 0.1, alpha = 1) +
  geom_errorbarh(data = sd_regions, aes(y = CAP2_mean, xmin = CAP1_mean - CAP1_sd, xmax = CAP1_mean + CAP1_sd, color = FarmNative_Crop), height = 0.1, alpha = 1) +
  # Add centroids
  geom_point(data = centroids, aes(x = CAP1_mean, y = CAP2_mean, color = FarmNative_Crop), size = 5, shape = 21, fill = "white", stroke = 1.5) +
  # Add labels for centroids
  geom_text(data = centroids, aes(x = CAP1_mean, y = CAP2_mean, label = FarmNative_Crop, color = FarmNative_Crop), 
            vjust = -1.2, hjust = 1.2, fontface = "bold", size = 4) +
  # Add significant environmental vectors
  geom_segment(data = biplot_scores, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "black", size = 0.7) +
  # Add labels for environmental vectors
  geom_text(data = biplot_scores, aes(x = CAP1, y = CAP2, label = Variable), 
            hjust = 1.2, vjust = 1.2, size = 4, color = "black") +
  scale_color_manual(values = group_colors) +
  labs(title = "CAP Ordination with Centroids, SD, and Envfit Vectors", 
       x = "CAP Axis 1", y = "CAP Axis 2") +
  theme_minimal()

# Save and display the plot
ggsave("/CAP_Centroids_SD_Envfit.pdf", plot = cap_plot, width = 10, height = 8,scale = 0.5)
print(cap_plot)

# -------------------------------------------------
# Step 8: Beta Dispersion Analysis - Separate Variables
# -------------------------------------------------

# Perform beta dispersion analysis for FarmNative_Crop
cat("Beta Dispersion for FarmNative_Crop:\n")
beta_disp_farmnative <- betadisper(bray_dist, sample_data$FarmNative_Crop)
beta_disp_test_farmnative <- permutest(beta_disp_farmnative, permutations = 999)
cat("P-Value for FarmNative_Crop:", beta_disp_test_farmnative$tab[1, "Pr(>F)"], "\n\n")

# Perform beta dispersion analysis for MAT
cat("Beta Dispersion for MAT:\n")
beta_disp_mat <- betadisper(bray_dist, sample_data$MAT)
beta_disp_test_mat <- permutest(beta_disp_mat, permutations = 999)
cat("P-Value for MAT:", beta_disp_test_mat$tab[1, "Pr(>F)"], "\n\n")

# Perform beta dispersion analysis for Site_ID
cat("Beta Dispersion for Site_ID:\n")
beta_disp_siteid <- betadisper(bray_dist, sample_data$Site_ID)
beta_disp_test_siteid <- permutest(beta_disp_siteid, permutations = 999)
cat("P-Value for Site_ID:", beta_disp_test_siteid$tab[1, "Pr(>F)"], "\n\n")

# Bin PC1 into quartiles
sample_data$PC1_binned <- cut(
  sample_data$PC1,
  breaks = quantile(sample_data$PC1, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)

# Ensure PC1_binned is a factor
sample_data$PC1_binned <- factor(sample_data$PC1_binned)

# Perform beta dispersion analysis for binned PC1
cat("Beta Dispersion for Binned PC1:\n")
beta_disp_pc1_binned <- betadisper(bray_dist, sample_data$PC1_binned)
beta_disp_test_pc1_binned <- permutest(beta_disp_pc1_binned, permutations = 999)
cat("P-Value for Binned PC1:", beta_disp_test_pc1_binned$tab[1, "Pr(>F)"], "\n")

