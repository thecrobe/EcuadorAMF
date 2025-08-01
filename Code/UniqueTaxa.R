library(igraph)
library(microViz)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)

##### Identify taxa exclusive to each landuse and crop  #####
load(file = "Data/Phyloseq/AMF_Ecuador_SSU.Rdata")
amf@sam_data$FarmNative_Crop<- paste(amf@sam_data$FarmNative, amf@sam_data$Crop)
amf@sam_data$FarmNative_Crop <- factor(amf@sam_data$FarmNative_Crop, levels=c("Farm Maize","Native Maize","Farm Potato","Native Potato"))
amf<-tax_fix(amf)


##### Exclusive ASVs for each crop in farms only 
afm_m<-psmelt(amf)
# Psmelt and keep only the farm samples
amf_m <- psmelt(amf) %>% filter(FarmNative == "Farm")
amf_potato <- amf_m %>%
  filter(Crop == "Potato" & Abundance != 0)
amf_maize <- amf_m %>%
  filter(Crop == "Maize" & Abundance != 0)

# Identify unique OTUs for 'Farm' compared to 'Native'
unique_otus_potato <- setdiff(unique(amf_potato$OTU), unique(amf_maize$OTU))

# Identify unique OTUs for 'Native' compared to 'Farm'
unique_otus_maize <- setdiff(unique(amf_maize$OTU), unique(amf_potato$OTU))

# Combine the exlusive OTUs into a data frame
unique_otus_counts <- data.frame(
  Category = c("Potato", "Maize"),
  Count = c(length(unique_otus_potato), length(unique_otus_maize))
)
