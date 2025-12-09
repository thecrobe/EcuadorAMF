library(phyloseq)
library(dplyr)
library(ggplot2)
library(microViz)

# === Load phyloseq object ===
load("Data/Phyloseq/AMF_Ecuador_SSU.Rdata")

# === Fix taxonomy ===
amf <- tax_fix(amf)
amf <- amf %>% tax_fix(unknowns = c(
  "Archaeospora Genus","Archaeosporaceae Family","Archaeosporales Order",
  "Claroideoglomus Genus","Glomeromycetes Class","Glomeromycotina sp.",
  "Mucoromycota Phylum","Mucoromycotina sp.","uncultured fungus",
  "uncultured Glomeromycotina","Ambispora Genus"
))

# === Create group ===
sample_data(amf)$FarmNative_Crop <- factor(
  paste(amf@sam_data$FarmNative, amf@sam_data$Crop),
  levels = c("Farm Maize", "Farm Potato", "Native Potato", "Native Maize")
)

ggsave(
  filename = "AMF_Relative_Abundance_Barplot.pdf",
  plot = p,
  width = 8,
  height = 6
)

df %>%
  arrange(FarmNative_Crop, desc(RelAbundance)) %>%
  mutate(RelAbundance = round(RelAbundance * 100, 2)) %>% print(n=1000)


df <- psmelt(amf) %>%
  mutate(
    Family = ifelse(
      Family %in% c("Acaulosporaceae","Ambisporaceae","Archaeosporaceae",
                    "Claroideoglomeraceae","Diversisporaceae","Gigasporaceae",
                    "Glomeraceae","Pacisporaceae","Paraglomeraceae"),
      Family,
      "Other"
    )
  ) %>%
  group_by(FarmNative_Crop, Family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(FarmNative_Crop) %>%
  mutate(RelAbundance = Abundance / sum(Abundance))


p <- ggplot(df, aes(x = FarmNative_Crop, y = RelAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "", y = "Relative abundance (%)", fill = "Family") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
