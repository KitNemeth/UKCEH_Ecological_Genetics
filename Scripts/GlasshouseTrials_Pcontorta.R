### Pinus contorta

GlasshouseTrial_Pcontorta <- read.delim("C:/Users/krinem/OneDrive - UKCEH/GlasshouseTrial_Pcontorta.txt")

library(ggpubr)
library(ggplot2)

#my_comparisons <- list( c("L", "M"), c("M", "H"), c("L", "H") )

GlasshouseTrial_Pcontorta$Adjusted.Height..mm. <- as.numeric(GlasshouseTrial_Pcontorta$Adjusted.Height..mm.)

GlasshouseTrial_Pcontorta$DBB..mm. <- as.numeric(GlasshouseTrial_Pcontorta$DBB..mm.)

my_comparisons <- list( c("Origin", "Plantation"), c("Origin", "Regeneration"), c("Plantation", "Regeneration") )


###Height 

OriginPlot_Height <- ggplot(GlasshouseTrial_Pcontorta, aes(x=Type, y=GlasshouseTrial_Pcontorta$Adjusted.Height..mm., fill=Type)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Seed source") + ylab("Height (mm)") +
  stat_compare_means(method= "t.test", comparisons = my_comparisons, label.y = c(80, 85, 90), hide.ns = TRUE) +
  stat_compare_means(method= "anova", label.y = 95)
OriginPlot_Height

##DBB
OriginPlot_DBB <- ggplot(GlasshouseTrial_Pcontorta, aes(x=Type, y=GlasshouseTrial_Pcontorta$DBB..mm., fill=Type)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Seed source") + ylab("DBB (mm)") +
  stat_compare_means(method= "t.test", comparisons = my_comparisons, label.y = c(0.9, 0.95, 1), hide.ns = TRUE) +
  stat_compare_means(method= "anova", label.y = 1.05)

OriginPlot_DBB

###Height 

OriginProvenancePlot_Height <- ggplot(GlasshouseTrial_Pcontorta, aes(x=Provenance, y=GlasshouseTrial_Pcontorta$Adjusted.Height..mm., fill=Type)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_brewer(palette="Set1") +
  xlab("Provenance") + ylab("Height (mm)") +
  stat_compare_means(method= "anova", label.y = 90) 

OriginProvenancePlot_Height

##DBB
OriginProvenancePlot_DBB <- ggplot(GlasshouseTrial_Pcontorta, aes(x=Provenance, y=GlasshouseTrial_Pcontorta$DBB..mm., fill=Type)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_brewer(palette="Set1") +
  xlab("Provenance") + ylab("DBB (mm)") +
  stat_compare_means(method= "anova", label.y = 1) 

OriginProvenancePlot_DBB

###Height 
my_comparisons <- list( c("Benmore", "Rowens"), c("Benmore", "NA"), c("Rowens", "NA") )

CollectionSitePlot_Height <- ggplot(GlasshouseTrial_Pcontorta_UKCollections, aes(x=GlasshouseTrial_Pcontorta_UKCollections$Collection.site..UK., y=GlasshouseTrial_Pcontorta_UKCollections$Adjusted.Height..mm., fill=GlasshouseTrial_Pcontorta_UKCollections$Collection.site..UK.)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Collection Site (UK)") + ylab("Height (mm)") +
  stat_compare_means(method = "t.test",label.x = 1.5, label.y = 80) 
CollectionSitePlot_Height

##DBB

CollectionSitePlot_DBB <- ggplot(GlasshouseTrial_Pcontorta_UKCollections, aes(x=GlasshouseTrial_Pcontorta_UKCollections$Collection.site..UK., y=GlasshouseTrial_Pcontorta_UKCollections$DBB..mm., fill=GlasshouseTrial_Pcontorta_UKCollections$Collection.site..UK.)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Collection Site (UK)") + ylab("DBB (mm)") +
  stat_compare_means(method = "t.test",label.x = 1.5, label.y = 1) 
CollectionSitePlot_DBB
