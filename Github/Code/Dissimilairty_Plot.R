# Load necessary libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggsignif)

# Extract OTU table and sample data
load(file = "AMF_Ecuador_SSU.Rdata")
otu_table_df <- as.data.frame(t(otu_table(amf)))
sample_data_df <- as.data.frame(sample_data(amf))

# Function to calculate Sørensen dissimilarity values
calculate_sorensen_values <- function(otu_table, sample_data, filter_condition) {
  # Filter sample data based on the provided condition
  samples_filtered <- sample_data[filter_condition, ]
  
  # Subset the OTU table to include only the filtered samples
  otu_table_filtered <- otu_table[rownames(samples_filtered), ]
  
  # Convert to presence/absence data
  otu_table_binary <- (otu_table_filtered > 0) * 1
  
  # Calculate Sørensen dissimilarity using betadiver
  dissimilarity_matrix <- betadiver(otu_table_binary, method = "sor")
  
  # Get the upper triangle of the matrix without the diagonal
  dissimilarity_values <- as.vector(dissimilarity_matrix[upper.tri(dissimilarity_matrix, diag = FALSE)])
  
  return(dissimilarity_values)
}

# Define filter conditions for Potato farms, Maize farms, and Uncultivated sites
filter_potato <- sample_data_df$Crop == "Potato" & sample_data_df$FarmNative == "Farm"
filter_maize <- sample_data_df$Crop == "Maize" & sample_data_df$FarmNative == "Farm"
filter_uncultivated <- sample_data_df$FarmNative == "Native"

# Calculate Sørensen dissimilarity values for Potato farms, Maize farms, and Uncultivated sites
sorensen_values_potato <- calculate_sorensen_values(otu_table_df, sample_data_df, filter_potato)
sorensen_values_maize <- calculate_sorensen_values(otu_table_df, sample_data_df, filter_maize)
sorensen_values_uncultivated <- calculate_sorensen_values(otu_table_df, sample_data_df, filter_uncultivated)

# Create a data frame for plotting
dissimilarity_df <- data.frame(
  Group = rep(c("Potato", "Maize", "Uncultivated"), 
              times = c(length(sorensen_values_potato), length(sorensen_values_maize), length(sorensen_values_uncultivated))),
  Dissimilarity = c(sorensen_values_potato, sorensen_values_maize, sorensen_values_uncultivated)
)

# Perform pairwise t-tests
pairwise_t_test_results <- pairwise.t.test(dissimilarity_df$Dissimilarity, dissimilarity_df$Group, p.adjust.method = "bonferroni")

# Extract p-values for the significant comparisons
pairwise_p_values <- as.data.frame(pairwise_t_test_results$p.value)

# Create a list of comparisons
comparisons <- list(
  c("Potato", "Maize"),
  c("Potato", "Uncultivated"),
  c("Maize", "Uncultivated")
)

# Boxplot of Sørensen dissimilarity values with jittered points
dissimplot <- ggplot(dissimilarity_df, aes(x = Group, y = Dissimilarity, color = Group)) +
  geom_jitter(color = "grey1", size = 0.5, width = 0.2, alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha = 1, lwd = 1, fill = NA) +
  scale_color_manual(values = c("Potato" = "#EA3929", "Maize" = "#EA3929", "Uncultivated" = "#1AB7D2")) +
  labs(title = "",
       x = "Land use",
       y = "Sørensen Dissimilarity") +
  theme_classic() +
  ylim(0, 0.7) +  # Adjust the y-axis limit to show from 0 to 0.7
  geom_signif(comparisons = comparisons, colour = "black",
              map_signif_level = TRUE, 
              test = "t.test", 
              y_position = c(0.52, 0.56, 0.6),  # Adjust these values based on your data
              tip_length = 0.01) 

# Save the plot to a file
##ggsave(filename = "Avg_Dissimilarity.pdf", plot = dissimplot)

