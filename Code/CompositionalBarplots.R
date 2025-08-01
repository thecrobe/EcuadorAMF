library(igraph)
library(microViz)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)

##### Identify changes in community composition
load(file = "Data/Phyloseq/AMF_Ecuador_SSU.Rdata")
amf@sam_data$FarmNative_Crop<- paste(amf@sam_data$FarmNative, amf@sam_data$Crop)
amf@sam_data$FarmNative_Crop <- factor(amf@sam_data$FarmNative_Crop, levels=c("Farm Maize","Native Maize","Farm Potato","Native Potato"))
amf<-tax_fix(amf)
amf1<-amf %>% tax_fix(unknowns = c("Archaeospora Genus", "Archaeosporaceae Family", "Archaeosporales Order", "Claroideoglomus Genus", "Glomeromycetes Class", "Glomeromycotina sp.", "Mucoromycota Phylum", "Mucoromycotina sp.", "uncultured fungus", "uncultured Glomeromycotina","Ambispora Genus", "Archaeospora Genus", "Archaeosporaceae Family", "Archaeosporales Order", "Claroideoglomus Genus", "Glomeromycetes Class", "Mucoromycota Phylum"))
library(phyloseq)
library(dplyr)
library(ggplot2)

# === Load phyloseq object ===
load("Data/Phyloseq/AMF_Ecuador_SSU.Rdata")

# === Create group ===
sample_data(amf)$FarmNative_Crop <- factor(
  paste(amf@sam_data$FarmNative, amf@sam_data$Crop),
  levels = c("Farm Maize", "Native Maize", "Farm Potato", "Native Potato")
)


# === Melt and process ===
df <- psmelt(amf) %>%
  mutate(Family = ifelse(is.na(Family) | Family == "", "Unknown", Family)) %>%
  group_by(FarmNative_Crop, Family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(FarmNative_Crop) %>%
  mutate(RelAbundance = Abundance / sum(Abundance))

# === Plot ===
p<-ggplot(df, aes(x = FarmNative_Crop, y = RelAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Group", y = "Relative Abundance", fill = "Family") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# === Identify most abundant family per group ===
top_families <- df %>%
  group_by(FarmNative_Crop) %>%
  slice_max(order_by = RelAbundance, n = 1) %>%
  arrange(FarmNative_Crop)

print(top_families)


ggsave(
  filename = "AMF_Relative_Abundance_Barplot.pdf",
  plot = p,
  width = 8,
  height = 6
)
