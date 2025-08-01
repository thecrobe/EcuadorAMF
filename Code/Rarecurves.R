library(MicrobiotaProcess)
library(cowplot)
library(dplyr)
library(phyloseq)
# for reproducibly random number
set.seed(1024)
crop_merged = merge_samples(amf, "Crop")
meta<-data.frame(crop_merged@sam_data)
crop_merged@sam_data$Crop<-row.names(meta)
rareres <- get_rarecurve(obj=crop_merged, chunks=400)

curve_crop<-ggrarecurve(obj=rareres, indexNames=c("Chao1")) +theme_cowplot(12) + 
  xlab("Sampling units") +  
  ylab("Richness")  +
  scale_color_discrete(name = "Crop") 

set.seed(1024)
meta$FarmNative
land_merged = merge_samples(amf, "FarmNative")
meta<-data.frame(land_merged@sam_data)
land_merged@sam_data$Land<-row.names(meta)
rareres <- get_rarecurve(obj=land_merged, chunks=400)

curve_land<-ggrarecurve(obj=rareres, indexNames=c("Chao1")) +theme_cowplot(12) + 
  xlab("Sampling units") + 
  ylab("Richness")  +
  scale_color_discrete(name = "Landuse") 

#All samples 
rareres <- get_rarecurve(obj=amf, chunks=400)
#all<-ggrarecurve(obj=rareres, indexNames=c("Chao1")) +theme_cowplot(12) + 
  xlab("Sampling units") + 
  ylab("Richness")  #+
  #scale_color_discrete(name = "Landuse") 

curves<-ggarrange(plotlist = list(curve_land,curve_crop))
ggsave(filename = "Rarefraction_curves.pdf",device = "pdf",scale = 0.8)

          