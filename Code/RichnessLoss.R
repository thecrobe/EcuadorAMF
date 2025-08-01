# Load necessary libraries for data processing, visualization, and modeling
library(dplyr)      # Data manipulation
library(tidyr)      # Data tidying
library(ggplot2)    # Plotting
library(ggsci)      # Scientific color schemes for ggplot
library(lme4)       # Linear mixed-effects models
library(MuMIn)      # Model dplyr::selection
library(purrr)      # Functional programming
library(phyloseq)   # Microbiome data manipulation
library(iNEXT)      # Diversity estimation (Hill numbers)
library(tidyverse)  # Data manipulation and visualization
library(brms)       # Bayesian modeling
library(gridExtra)  # For arranging plots

# Set seed for reproducibility
set.seed(777)

# Step 1: Load and clean phyloseq data
setwd("/")
load("Data/Phyloseq/AMF_Ecuador_SSU.Rdata")  # Load phyloseq object
amf <- prune_taxa(taxa_sums(amf) > 0, amf)    # Remove taxa with zero counts

# Step 2: Calculate Hill Numbers (Species Richness)
iNEXT_input_table <- otu_table(amf) %>%
  data.frame()
iNEXT_data <- iNEXT(iNEXT_input_table, q = c(0,1), datatype = "abundance")
richness <- iNEXT_data$AsyEst %>%
  filter(Diversity == "Species richness")
amf@sam_data$Chao1 <- richness$Estimator  # Add richness estimates to sample data
mapping <- data.frame(amf@sam_data)

# Step 3: Calculate Mean Richness by Land Use and Crop Type
means <- mapping %>%
  group_by(FarmNative, Crop) %>%
  summarize(Mean_Chao1 = mean(round(Chao1), na.rm = TRUE), .groups = 'drop') %>%
  mutate(Mean_Chao1 = round(Mean_Chao1))
print(means)

# Step 4: Calculate Percent Change in Richness
farms <- mapping %>%
  filter(FarmNative == "Farm") %>%
  dplyr::select(Site_ID, Chao1_Farm = Chao1)
natives <- mapping %>%
  filter(FarmNative == "Native") %>%
  dplyr::select(Site_ID, Crop, Chao1_Native = Chao1)
comparison <- left_join(farms, natives, by = "Site_ID")

compare <- comparison %>%
  mutate(Percent_Change = ((Chao1_Native - Chao1_Farm) / Chao1_Farm) * 100)
summary(compare$Percent_Change)

# Step 5: Visualize Percent Change in Richness
compare_plot <- ggplot(compare, aes(x = Percent_Change, y = "", color = Crop)) +
  geom_jitter(height = 0.3, size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c("#CCCC00","#AE64CC")) +
  theme_minimal() +
  xlab("Richness (% difference)") +
  ylab("")

# Step 6: Add variables from remote sensing
data <- read.csv("Data/EcuadorLayers_GEE.csv") %>%
  left_join(mapping, .) %>%
  rename(MAP = CHELSA_BIO_Annual_Precipitation, MAT = CHELSA_BIO_Annual_Mean_Temperature) %>%
  mutate(MAT = MAT * 0.1, MAP = MAP * 0.1, IsFarm = as.factor(FarmNative == "Farm"))

# Step 7: Calculate PCA for Soil Chemistry Data
data$CN <- data$Carbon_totalOrganic_percent / data$Nitrogen_percent
data$NP <- data$Nitrogen_percent / data$P_mgL
chem <- data %>%
  dplyr::select(Crop, CN, NP, pH, FarmNative) %>%
  mutate(CN = CN * 100)
pca_result <- prcomp(chem %>% dplyr::select(where(is.numeric)), center = TRUE, scale. = TRUE)

# Step 8: Prepare Data for Bayesian Models
df_scaled <- data %>%
  dplyr::select(MAP, MAT, NP, CN, pH) %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  mutate(Site_ID = data$Site_ID, Crop = data$Crop, PC1 = pca_result$x[,1])
df_scaled$Chao1 <- data$Chao1
df_scaled$IsFarm <- data$IsFarm

# Step 9: Define models with and without `IsFarm` as a random slope
model_list <- list(
#  "Intercept_no_slope" = "Chao1 ~ 1 + (1 | Site_ID)",
#  "Soil_no_slope" = "Chao1 ~ CN + NP + pH + (1 | Site_ID)",
#  "Crop_no_slope" = "Chao1 ~ Crop + (1 | Site_ID)",
#  "Landuse_no_slope" = "Chao1 ~ IsFarm + (1 | Site_ID)",
#  "Climate_no_slope" = "Chao1 ~ MAP + MAT + (1 | Site_ID)",
#  "ClimateLanduse_add_no_slope" = "Chao1 ~ MAP + MAT + IsFarm + (1 | Site_ID)",
#  "ClimateLanduse_add_soil_no_slope" = "Chao1 ~ MAP + MAT + IsFarm + CN + NP + pH + (1 | Site_ID)",
#  "ClimateLanduse_int_no_slope" = "Chao1 ~ MAT * IsFarm + MAP * IsFarm + CN + NP + pH + (1 | Site_ID)",
  
#  "Intercept_with_slope" = "Chao1 ~ 1 + (1 + IsFarm | Site_ID)",
#  "Soil_with_slope" = "Chao1 ~ CN + NP + pH + (1 + IsFarm | Site_ID)",
#  "Crop_with_slope" = "Chao1 ~ Crop + (1 + IsFarm | Site_ID)",
#  "Landuse_with_slope" = "Chao1 ~ IsFarm + (1 + IsFarm | Site_ID)",
#  "Climate_with_slope" = "Chao1 ~ MAP + MAT + (1 + IsFarm | Site_ID)",
#  "ClimateLanduse_add_with_slope" = "Chao1 ~ MAP + MAT + IsFarm + (1 + IsFarm | Site_ID)",
#  "ClimateLanduse_add_soil_with_slope" = "Chao1 ~ MAP + MAT + IsFarm + CN + NP + pH + (1 + IsFarm | Site_ID)",
  "ClimateLanduse_int_with_slope" = "Chao1 ~ MAT * IsFarm + MAP * IsFarm + CN + NP + pH + (1 + IsFarm | Site_ID)"
)

# Step 10: Fit models and store WAIC and R² values
fit_and_evaluate <- function(formula, data) {
  model <- brm(formula, data = data, refresh = 0, sample_prior = "yes")
  waic_val <- waic(model)$estimates["waic", "Estimate"]
  r2_val <- bayes_R2(model)[1]
  return(list(waic = waic_val, r2 = r2_val))
}

results <- map_dfr(model_list, ~fit_and_evaluate(.x, df_scaled), .id = "Model")
results <- results %>%
  mutate(Slope_Type = if_else(grepl("no_slope", Model), "Without Random Slope", "With Random Slope"))

# Step 11: Calculate ΔWAIC from the intercept model
intercept_waic_no_slope <- results %>% filter(Model == "Intercept_no_slope") %>% pull(waic)
intercept_waic_with_slope <- results %>% filter(Model == "Intercept_with_slope") %>% pull(waic)

results <- results %>%
  mutate(Delta_WAIC = ifelse(Slope_Type == "Without Random Slope",
                             waic - intercept_waic_no_slope,
                             waic - intercept_waic_with_slope))

# Step 12: Plot ΔWAIC and R² values
delta_waic_plot <- ggplot(results, aes(x = Model, y = Delta_WAIC, fill = Slope_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "ΔWAIC Values for Models", x = "Model", y = "ΔWAIC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

r2_plot <- ggplot(results, aes(x = Model, y = r2, fill = Slope_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "R² Values for Models", x = "Model", y = "R²") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot <- ggpubr::ggarrange(delta_waic_plot, r2_plot, 
                                   ncol = 2, common.legend = TRUE, legend = "right")

print(combined_plot)
#ggsave(plot = combined_plot, filename = "", scale = 0.8)

# Step 13: Conditional effects and hypothesis testing
best_fit <- brm(
  Chao1 ~ MAT * IsFarm + MAP * IsFarm + CN + NP + pH + (1 + IsFarm | Site_ID),
  data = df_scaled,
  iter = 10000,
  warmup = 2000,
  chains = 4,
  control = list(adapt_delta = 0.99),
  sample_prior = TRUE
)


effects <- conditional_effects(best_fit, effects = "MAT:IsFarm", prob = 0.95)
effects_df <- as.data.frame(effects$`MAT:IsFarm`)
MATxFarm_Plot <- ggplot(effects_df, aes(x = MAT, y = estimate__, color = IsFarm)) +
  stat_smooth(method = "lm", se = FALSE, size = 3) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.1) +
  theme_minimal()

effects <- conditional_effects(best_fit, effects = "MAP:IsFarm", prob = 0.95)
effects_df <- as.data.frame(effects$`MAP:IsFarm`)
MAPxFarm_Plot <- ggplot(effects_df, aes(x = MAP, y = estimate__, color = IsFarm)) +
  stat_smooth(method = "lm", se = FALSE, size = 3) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.1) +
  theme_minimal()

multi <- ggpubr::ggarrange(MATxFarm_Plot, MAPxFarm_Plot, ncol = 2, legend = FALSE)
#ggsave(plot = multi, filename = "", scale = 0.8, width = 15, height = 6)

# Additional hypothesis testing and model summaries
# Standard deviations (if needed later for scaling)
sd_mat <- 2.58     # Standard deviation for MAT (degrees Celsius)
sd_map <- 17.49    # Standard deviation for MAP (mm precipitation)

# Extract coefficients and uncertainties
coef_summary <- posterior_summary(best_fit)

# MAT results for IsFarm = FALSE (uncultivated)
mat_coef_farm_false <- coef_summary["b_MAT", "Estimate"]
mat_uncertainty_farm_false <- coef_summary["b_MAT", "Est.Error"]
best_fit
# Print MAT results for IsFarm = FALSE
print("Coefficient for MAT when IsFarm = FALSE (uncultivated):")
print(mat_coef_farm_false)
print("Uncertainty (Est.Error) for MAT when IsFarm = FALSE (uncultivated):")
print(mat_uncertainty_farm_false)

# MAT hypothesis testing
hypothesis_false <- hypothesis(best_fit, "MAT > 0")
hypothesis_true <- hypothesis(best_fit, "MAT:IsFarmTRUE > 0")

# Print hypothesis testing results for MAT
print("Hypothesis test for MAT when IsFarm = FALSE (uncultivated):")
print(hypothesis_false)
print("Hypothesis test for MAT when IsFarm = TRUE (farm):")
print(hypothesis_true)
hypothesis_gradients <- hypothesis(best_fit, "MAT:IsFarmTRUE < 2.27",alpha = 0.05) #test if climate gradients sig,NA is significant


# MAP results for IsFarm = FALSE (uncultivated)
map_coef_farm_false <- coef_summary["b_MAP", "Estimate"]
map_uncertainty_farm_false <- coef_summary["b_MAP", "Est.Error"]

# Print MAP results for IsFarm = FALSE
print("Coefficient for MAP when IsFarm = FALSE (uncultivated):")
print(map_coef_farm_false)
print("Uncertainty (Est.Error) for MAP when IsFarm = FALSE (uncultivated):")
print(map_uncertainty_farm_false)

# MAP hypothesis testing
hypothesis_map_false <- hypothesis(best_fit, "MAP < 0")
hypothesis_map_true <- hypothesis(best_fit, "IsFarmTRUE:MAP > 0")

# Print hypothesis testing results for MAP
print("Hypothesis test for MAP when IsFarm = FALSE (uncultivated):")
print(hypothesis_map_false)
print("Hypothesis test for MAP when IsFarm = TRUE (farm):")
print(hypothesis_map_true)
hypothesis_gradients <- hypothesis(best_fit, "IsFarmTRUE:MAP = -4.71") #test if climate gradients sig differ


# Identify taxa exclusive to each land use and crop
library(igraph)
library(microViz)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)

# Step 14: Prepare data for land use analysis
load("Data/Phyloseq/AMF_Ecuador_SSU.Rdata")
amf@sam_data$FarmNative_Crop <- paste(amf@sam_data$FarmNative, amf@sam_data$Crop)
amf@sam_data$FarmNative_Crop <- factor(amf@sam_data$FarmNative_Crop, levels = c("Farm Maize", "Native Maize", "Farm Potato", "Native Potato"))
amf <- tax_fix(amf)

# Step 15: Calculate Hill Numbers using iNEXT and add to phyloseq object
iNEXT_input_table <- otu_table(amf) %>%
  data.frame()
iNEXT_data <- iNEXT(iNEXT_input_table, q = c(0, 1), datatype = "abundance")
richness <- iNEXT_data$AsyEst %>% filter(Diversity == "Species richness")
amf@sam_data$Chao1 <- richness$Estimator
mapping <- data.frame(amf@sam_data)

# Step 16: Identify exclusive taxa for each land use
amf_m <- psmelt(amf)
amf_farm <- amf_m %>% filter(FarmNative == "Farm" & Abundance != 0)
amf_native <- amf_m %>% filter(FarmNative == "Native" & Abundance != 0)

# Find unique OTUs in each land use category
unique_otus_farm <- setdiff(unique(amf_farm$OTU), unique(amf_native$OTU))
unique_otus_native <- setdiff(unique(amf_native$OTU), unique(amf_farm$OTU))

# Combine exclusive OTUs into a data frame for visualization
unique_otus_counts <- data.frame(
  Category = c("Monoculture", "Uncultivated"),
  Count = c(length(unique_otus_farm), length(unique_otus_native))
)

plot_unique_land <- ggplot(unique_otus_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(x = "Land Use", y = "Number of ASVs") +
  ggsci::scale_fill_npg() + 
  theme_minimal() +
  ylim(0, 400)

# Step 17: Create phyloseq object with only exclusive OTUs
merged_otus <- union(unique_otus_farm, unique_otus_native)
amf_subset <- prune_taxa(merged_otus, amf)

# List taxa to remove (optional, for improved clarity in plots)
taxa_to_remove <- c("Ambispora", "Archaeospora", "Archaeosporaceae Family", "Archaeosporales Order")
ps_filtered_exclusive <- subset_taxa(amf_subset, !Genus %in% taxa_to_remove)

# Define custom color palette for plotting
myPal <- tax_palette(data = ps_filtered_exclusive, rank = "Genus", n = 15, pal = "greenArmytage", add = c(Other = "white"))

# Step 18: Plot exclusive taxa composition by land use
plot_unique_taxa <- ps_filtered_exclusive %>%
  tax_fix(unknowns = c("Archaeospora Genus", "Archaeosporaceae Family", "Archaeosporales Order", "Glomeromycotina sp.", "uncultured Glomeromycotina", "uncultured Glomus")) %>%
  phyloseq::merge_samples(group = "FarmNative") %>%
  comp_barplot(tax_level = "Genus", palette = myPal, taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), 
               other_name = "Other taxa", n_taxa = 10, bar_width = 0.8) + 
  xlab("Land Use") + ggtitle("Land Use Exclusive Communities")

# Step 19: Run PERMANOVA to test for differences in unique taxa across land uses
ps_filtered_exclusive <- prune_taxa(taxa_sums(ps_filtered_exclusive) > 0, ps_filtered_exclusive)
otu_table <- data.frame(t(otu_table(ps_filtered_exclusive)))
sample_data <- data.frame(sample_data(ps_filtered_exclusive))
permanova_result <- adonis2(otu_table ~ FarmNative, data = sample_data, permutations = 999)
print(permanova_result)

# Step 20: Identify unique ASVs for each crop type in farm samples
amf_m <- psmelt(amf) %>% filter(FarmNative == "Farm")
amf_potato <- amf_m %>% filter(Crop == "Potato" & Abundance != 0)
amf_maize <- amf_m %>% filter(Crop == "Maize" & Abundance != 0)

unique_otus_potato <- setdiff(unique(amf_potato$OTU), unique(amf_maize$OTU))
unique_otus_maize <- setdiff(unique(amf_maize$OTU), unique(amf_potato$OTU))

# Combine exclusive OTUs by crop type into a data frame
unique_otus_counts_crop <- data.frame(
  Category = c("Potato", "Maize"),
  Count = c(length(unique_otus_potato), length(unique_otus_maize))
)

plot_unique_crop <- ggplot(unique_otus_counts_crop, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(x = "Crop", y = "Number of ASVs") +
  scale_fill_manual(values = c("#CCCC00", "#AE64CC")) +
  theme_minimal() +
  ylim(0, 400)

# Step 21: Plot unique taxa composition by crop
AMF_farm <- subset_samples(amf, FarmNative == "Farm")
merged_otus_farms <- union(unique_otus_potato, unique_otus_maize)
amf_subset <- prune_taxa(merged_otus_farms, AMF_farm)

ps_filtered_exclusive_crop <- subset_taxa(amf_subset, !Genus %in% taxa_to_remove)
plot_unique_taxa_crop <- ps_filtered_exclusive_crop %>%
  tax_fix(unknowns = c("Archaeospora Genus", "Archaeosporaceae Family", "Archaeosporales Order", "Glomeromycotina sp.", "uncultured Glomeromycotina", "uncultured Glomus")) %>%
  phyloseq::merge_samples(group = "Crop") %>%
  comp_barplot(tax_level = "Genus", palette = myPal, taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), 
               other_name = "Other taxa", n_taxa = 10, bar_width = 0.8) + 
  xlab("Crop") + ggtitle("Crop Exclusive Communities")

# Step 22: Run PERMANOVA to test for differences in unique taxa across crops
otu_table <- data.frame(t(otu_table(ps_filtered_exclusive_crop)))
sample_data <- data.frame(sample_data(ps_filtered_exclusive_crop))
permanova_result_crop <- adonis2(otu_table ~ Crop, data = sample_data, permutations = 999)
print(permanova_result_crop)

# Step 23: Create multi-panel plot for land use and crop exclusive taxa
multi_panel_plot <- ggpubr::ggarrange(
  plot_unique_land, plot_unique_crop, plot_unique_taxa, plot_unique_taxa_crop,
  labels = "auto", ncol = 2, nrow = 2, common.legend = TRUE
)

# Save the multi-panel plot
#ggsave(plot = multi_panel_plot, filename = "", scale = 0.8, width = 15, height = 10)
