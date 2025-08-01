######### This script is for plotting and describing the sample distributions per ecoregion
######### Author: Justin D Stewart
######### Date: Dec 24, 2022

#################### Readins
###Libraries
library(tidyverse)
library(raster)
library(cowplot)
library(spdplyr)
library(broom)
library(reshape2)

library(tidyverse)
library(sf)
library(sp)
library(rnaturalearth)
library(rmapshaper)
library(terra)
library(raster)
library(viridis)

#load("Data/Phyloseq/AMF_ecuador.rds") #load in phyloseq object
mapping<-data.frame(sample_data(amf))
potato<-mapping %>% filter(Crop=="Potato")
maize<-mapping %>% filter(Crop=="Maize")



# Plot the map
map <- ggplot() +
  geom_sf(data = ecuador, fill = "gray90", color = "gray10") +
  geom_sf(data = potato_sf, aes(color = "Potato"), size = 3, shape = 21, fill = "#915AA2") +
  geom_sf(data = maize_sf, aes(color = "Maize"), size = 3, shape = 21, fill = "#C1C72B") +
  scale_color_manual(values = c("Potato" = "#915AA2", "Maize" = "#C1C72B")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray80"),
    panel.background = element_rect(fill = "white")
  ) +
  labs(title = "",
       x = "Longitude",
       y = "Latitude")

##### CLimate plots 
# Read in additional environmental data
layers <- read.csv("Data/EcuadorLayers_GEE.csv", header = TRUE)

# Assuming 'amf' and 'sample_data' are already loaded
sample_data <- data.frame(sample_data(amf))
sample_data <- left_join(sample_data, layers, by = "Site_ID")

# Rename and adjust environmental variables
sample_data <- sample_data %>% 
  rename(MAP = CHELSA_BIO_Annual_Precipitation,
         MAT = CHELSA_BIO_Annual_Mean_Temperature) %>%
  mutate(MAT = MAT * 0.1,
         MAP = MAP * 0.1,
         P_percent = (P_mgL / 1000) * 10,
         CN = Carbon_totalOrganic_percent / (Nitrogen_percent * 10000),
         NP = (Nitrogen_percent * 10000) / P_percent)

# Filter data for Potato and Maize
potato_data <- sample_data %>% filter(Crop == "Potato")
maize_data <- sample_data %>% filter(Crop == "Maize")

# Define color for crops
crop_colors <- c("Potato" = "#915AA2", "Maize" = "#C1C72B")

# Function to plot density plots with specified x-axis range and labels, rescaled to 0-1
plot_density <- function(data1, data2, variable, x_min, x_max, x_label) {
  ggplot() +
    geom_density(data = data1, aes_string(x = variable, color = "Crop", y = "..scaled.."), size = 1.5) +
    geom_density(data = data2, aes_string(x = variable, color = "Crop", y = "..scaled.."), size = 1.5) +
    scale_color_manual(values = crop_colors) +
    theme_minimal() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray10"),
      panel.background = element_rect(fill = "gray95")
    ) +
    labs(title = "", x = x_label, y = "Density") +
    xlim(x_min, x_max) +
    ylab("") + theme_minimal()
}

# Create density plots with specified x-axis ranges and labels
plot1 <- plot_density(potato_data, maize_data, "MAP", 50, 450, expression(MAP ~ (kg ~ m^{-2})))
plot2 <- plot_density(potato_data, maize_data, "MAT", 5, 20, expression(MAT ~ (degree ~ C)))

# Arrange the plots into a single column
combined_plot <- ggarrange(plot1, plot2, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

# Combine the map and density plots using patchwork
final_plot <- map + (plot1 / plot2) + plot_layout(ncol = 2, widths = c(2, 1))

# Display the final combined plot and save 
print(final_plot)
output_path <- "/Fig1_map.pdf"
ggsave(output_path, plot = final_plot, device = "pdf", width = 10, height = 7)
