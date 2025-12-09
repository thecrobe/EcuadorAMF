# Load necessary libraries
library(microViz)
library(emmeans)
library(reshape2)
library(tidyverse)
library(sp)
library(spdep)
library(phyloseq)
library(ggpubr)
library(vegan)
library(betareg)
library(modelbased)
library(see)
library(dplyr)
library(brms)
library(ggsci)

##### Step 1: Data Preparation #####

# Set working directory and load phyloseq object
load("Data/Phyloseq/AMF_Ecuador_SSU.Rdata")


# Clean and format phyloseq data
amf <- tax_fix(amf)
amf <- phyloseq_validate(amf, remove_undetected = TRUE)
mapping <- data.frame(sample_data(amf))
mapping <- mapping[!duplicated(colnames(mapping))]

# Import remote sensing layers, join, rename, and scale values
layers <- read.csv("Data/EcuadorGEE.csv", header = TRUE)
data <- left_join(mapping, layers) %>%
  rename(MAP = CHELSA_BIO_Annual_Precipitation, MAT = CHELSA_BIO_Annual_Mean_Temperature) %>%
  mutate(
    MAT = MAT * 0.1,
    MAP = MAP * 0.1,
    P_percent = (P_mgL / 1000) * 10,
    CN = Carbon_totalOrganic_percent / (Nitrogen_percent * 10000),
    NP = (Nitrogen_percent * 10000) / P_percent,
    IsFarm = as.factor(ifelse(FarmNative == "Farm", 1, 0)),
    IsMaize = as.factor(ifelse(Crop == "Maize", 1, 0))
  )

# Summarize soil chemistry per group
chem <- data %>%
  dplyr::select(Crop, CN, NP, pH, FarmNative) %>%
  mutate(CN = CN * 100)
result <- chem %>%
  group_by(FarmNative, Crop) %>%
  summarize(across(c(CN, NP, pH), list(mean = mean, sd = sd), .names = "{col}_{fn}"))

##### Step 2: Calculate Paired Dissimilarity #####

# Prepare OTU table for beta diversity calculation
mat <- t(data.frame(otu_table(amf)))
mat <- data.matrix(mat)
mat[mat > 0] <- 1

# Calculate replacement and nestedness components using betapart
library(betapart)
matrices <- beta.pair(mat, index.family = "sorensen")
replacement <- matrices$beta.sim
nestedness <- matrices$beta.sne

# Melt replacement and nestedness matrices, then filter allowed pairs
melted_replacement <- reshape2::melt(data.matrix(replacement))
melted_nestedness <- reshape2::melt(data.matrix(nestedness))

allowed_pairs <- tibble(
  Var1 = c("M1F", "M2F", "M3F", "M7F", "M8F", "M9F", "M10F", "M15F", "M18F", "M19F",
           "P11F", "P12F", "P13F", "P14F", "P15F", "P16F", "P17F", "P18F", "P19F",
           "P1F", "P3F", "P4F", "P5F", "P6F", "P7F", "P8F", "P9F"),
  Var2 = c("M1N1", "M2N1", "M3N1", "M7N1", "M8N1", "M9N1", "M10N1", "M15N1", "M18N1", "M19N1",
           "P11N1", "P12N1", "P13N1", "P14N1", "P15N1", "P16N1", "P17N1", "P18N1", "P19N1",
           "P1N1", "P3N1", "P4N1", "P5N1", "P6N1", "P7N1", "P8N1", "P9N1")
)

# Filter pairs and join replacement and nestedness matrices into paired_dissimilarity
paired_replacement <- melted_replacement %>%
  semi_join(allowed_pairs, by = c("Var1", "Var2")) %>%
  mutate(Site_ID = str_remove(Var1, "F"), replacement = value) %>%
  dplyr::select(Site_ID, replacement)

paired_nestedness <- melted_nestedness %>%
  semi_join(allowed_pairs, by = c("Var1", "Var2")) %>%
  mutate(Site_ID = str_remove(Var1, "F"), nestedness = value) %>%
  dplyr::select(Site_ID, nestedness)

paired_dissimilarity <- inner_join(paired_replacement, paired_nestedness, by = "Site_ID")

# Create soil ratio data and join with paired dissimilarity
variables_to_calculate_ratio <- c("CN", "NP", "pH")
soil_ratio <- data %>%
  group_by(Site_ID) %>%
  summarize(across(all_of(variables_to_calculate_ratio), list(Farm_to_Native_LogEffect = ~ log(.[FarmNative == "Farm"] / .[FarmNative == "Native"]))))

paired_beta <- inner_join(paired_dissimilarity, soil_ratio, by = "Site_ID")

# Calculate initial residual as the shared portion of the community
paired_beta$residual <- 1 - (paired_beta$nestedness + paired_beta$replacement)

# Define a small adjustment to prevent exact 0 or 1 values
small_adjustment <- 1e-10

# Apply adjustments to ensure no component is exactly 0 or 1, and all sum to 1
paired_beta <- paired_beta %>%
  rowwise() %>%
  mutate(
    # Ensure no component is zero by adding a small adjustment
    nestedness = nestedness + small_adjustment,
    replacement = replacement + small_adjustment,
    residual = residual + small_adjustment,
    
    # Recalculate total with adjustments
    total_sum = nestedness + replacement + residual,
    
    # Normalize components to ensure they sum exactly to 1
    nestedness = nestedness / total_sum,
    replacement = replacement / total_sum,
    residual = residual / total_sum,
    
    # Final check: if sum is slightly off due to rounding, adjust residual
    residual = 1 - (nestedness + replacement)
  ) %>%
  ungroup() %>%
  dplyr::select(-total_sum)  # Remove the temporary column

paired_beta <- paired_beta %>%
  left_join(data %>% dplyr::select(Site_ID, MAT, MAP), by = "Site_ID")


##### Step 3: Define and Fit Models #####

# Fit four models: intercept-only, temperature-only, precipitation-only, and temperature + precipitation
models <- list(
  "Intercept-only" = brm(cbind(residual, replacement, nestedness) ~ 1, data = paired_beta, 
                         family = dirichlet(), chains = 4, iter = 3000, warmup = 500, thin = 4, save_pars = save_pars(all = TRUE)),
  
  "Temperature-only (MAT)" = brm(cbind(residual, replacement, nestedness) ~ MAT, data = paired_beta, 
                                 family = dirichlet(), chains = 4, iter = 3000, warmup = 500, thin = 4, save_pars = save_pars(all = TRUE)),
  
  "Precipitation-only (MAP)" = brm(cbind(residual, replacement, nestedness) ~ MAP, data = paired_beta, 
                                   family = dirichlet(), chains = 4, iter = 3000, warmup = 500, thin = 4, save_pars = save_pars(all = TRUE)),
  
  "Temperature + Precipitation" = brm(cbind(residual, replacement, nestedness) ~ MAT + MAP, data = paired_beta, 
                                      family = dirichlet(), chains = 4, iter = 3000, warmup = 500, thin = 4, save_pars = save_pars(all = TRUE))
)


# ===============================================================
# Extract posterior summaries (fixed effects) for each model
# ===============================================================

coef_tables <- lapply(models, function(mod) {
  posterior_summary(mod, pars = "^b_") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter")
})

# Add model names
names(coef_tables) <- names(models)

# Combine into one data frame
coef_df <- dplyr::bind_rows(
  lapply(names(coef_tables), function(nm) {
    coef_tables[[nm]] %>% mutate(Model = nm)
  })
)

# Save
write.csv(coef_df, "//Users/justinstewart/Dropbox/Collaborations/Ecuador/Writing/MolEcol/Rev1/DirichletModel_Estimates.csv", row.names = FALSE)

# Preview
head(coef_df)


# Calculate WAIC and LOO for each model
waic_values <- sapply(models, function(model) waic(model)$estimates["waic", "Estimate"])
loo_values <- lapply(models, loo)

# Extract the WAIC of the intercept-only (null) model
null_waic <- waic_values["Intercept-only"]

# Calculate WAIC differences from the null model
waic_diff <- waic_values - null_waic

# Create a data frame for WAIC differences
waic_df <- data.frame(Model = names(waic_diff), WAIC_Diff = waic_diff) %>%
  arrange(desc(WAIC_Diff))  # Arrange models by WAIC difference from high to low

# LOO model comparison (optional for additional comparison)
loo_comparison <- loo_compare(loo_values)

##### Plot WAIC Differences from the Null Model #####

waic_plot <- ggplot(waic_df, aes(y = reorder(Model, WAIC_Diff), x = WAIC_Diff, fill = Model)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Model Comparison by WAIC Difference from Intercept model", y = "Model", x = "WAIC Difference") +
  theme_minimal() +
  theme(legend.position = "none")

# Display the plot
print(waic_plot)
#ggsave(filename = "DirichletModel_WAIC.png",scale = 0.8)

##### Step 5: Plot Marginal Effects #####

# Plot marginal effects for MAP and MAT from the Temperature-only and Precipitation-only models

# Marginal effects for MAT
MAT_effect <- marginal_effects(models[["Temperature + Precipitation"]], effects = "MAT", categorical = TRUE)
df_mat <- as.data.frame(MAT_effect$`MAT:cats__`)
mat_plot <- ggplot(df_mat, aes(x = MAT, y = estimate__, fill = cats__)) +
  geom_area() +
  ylim(0, 1) +
  xlab(expression("Temperature ("*degree*"C)")) +
  ylab("Sorenson dissimilarity") +
  scale_fill_manual(values = c("#1E88E5", "#FFC107", "grey40"), name = "Component") +
  theme_minimal()

# Marginal effects for MAP
MAP_effect <- marginal_effects(models[["Temperature + Precipitation"]], effects = "MAP", categorical = TRUE)
df_map <- as.data.frame(MAP_effect$`MAP:cats__`)
map_plot <- ggplot(df_map, aes(x = MAP, y = estimate__, fill = cats__)) +
  geom_area() +
  ylim(0, 1) +
  ylab("Sorenson dissimilarity") +
  xlab(expression("MAP (mm/m"^2*")")) +
  scale_fill_manual(values = c("#1E88E5", "#FFC107", "grey40"), name = "Component") +
  theme_minimal()

# Arrange MAT and MAP plots side by side
combined_plot <- ggarrange(mat_plot, map_plot, ncol = 2, labels = c("Temperature effect", "Precipitation effect"),common.legend = TRUE)

# Display the combined plot
print(combined_plot)

#ggsave("DirichletModel_plot.pdf",device = "pdf",plot = combined_plot,scale = 0.6,width = 15,height = 7)

#Get average effects for reporting. 
# Average effect for MAP
MAP_effects <- conditional_effects(models[["Temperature + Precipitation"]], effects = "MAP",categorical = TRUE)
MAP_avg_effect <- as.data.frame(MAP_effects$`MAP:cats__`)
summary(MAP_avg_effect$nestedness)
summary(MAP_avg_effect$replacement)
summary(MAP_avg_effect$residual)

# Average effect for MAT
mat_effects <- conditional_effects(models[["Temperature + Precipitation"]], effects = "MAT",categorical = TRUE)
mat_avg_effect <- as.data.frame(mat_effects$`MAT:cats__`)
summary(mat_avg_effect$nestedness)
summary(mat_avg_effect$replacement)
summary(mat_avg_effect$residual)


# Hypothesis testing
# Preform hypothesis testing to estimate probability of a significant effect. 
mod_chemistry_climate_add<-models$`Temperature + Precipitation`
best_fit_df<-data.frame(mod_chemistry_climate_add)
hypothesis(mod_chemistry_climate_add ,hypothesis = "munestedness_MAT > 0")
hypothesis(mod_chemistry_climate_add ,hypothesis = "mureplacement_MAT > 0")
hypothesis(mod_chemistry_climate_add ,hypothesis = "munestedness_MAP > 0")
hypothesis(mod_chemistry_climate_add ,hypothesis = "mureplacement_MAP > 0")
