### Tsuga heterophylla

GlasshouseTrial_Theterophylla <- read.delim("C:/Users/krinem/OneDrive - UKCEH/GlasshouseTrial_Theterophylla.txt")

library(ggpubr)
library(ggplot2)

GlasshouseTrial_Theterophylla$Height..mm. <- as.numeric(GlasshouseTrial_Theterophylla$Height..mm.)

GlasshouseTrial_Theterophylla$DBB..mm. <- as.numeric(GlasshouseTrial_Theterophylla$DBB..mm.)

GlasshouseTrial_TheterophyllaExtra <- GlasshouseTrial_Theterophylla[GlasshouseTrial_Theterophylla$Trial == "Extra", ]

GlasshouseTrial_Theterophylla <- GlasshouseTrial_Theterophylla[GlasshouseTrial_Theterophylla$Trial == "Main", ]


###Height 

OriginPlot_Height <- ggplot(GlasshouseTrial_Theterophylla, aes(x=Type, y=GlasshouseTrial_Theterophylla$Height..mm., fill=Type)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Seed source") + ylab("Height (mm)") +
  stat_compare_means(method = "t.test",label.x = 1.5, label.y = 45) 

OriginPlot_Height

##DBB
OriginPlot_DBB <- ggplot(GlasshouseTrial_Theterophylla, aes(x=Type, y=GlasshouseTrial_Theterophylla$DBB..mm., fill=Type)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Seed source") + ylab("DBB (mm)") +
  stat_compare_means(method = "t.test",label.x = 1.5, label.y = 1) 

OriginPlot_DBB
###Height 

OriginProvenancePlot_Height <- ggplot(GlasshouseTrial_Theterophylla, aes(x=Provenance, y=GlasshouseTrial_Theterophylla$Height..mm., fill=Type)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_brewer(palette="Set1") +
  xlab("Provenance") + ylab("Height (mm)") +
  stat_compare_means(method= "anova", label.y = 60)
OriginProvenancePlot_Height

##DBB
OriginProvenancePlot_DBB <- ggplot(GlasshouseTrial_Theterophylla, aes(x=Provenance, y=GlasshouseTrial_Theterophylla$DBB..mm., fill=Type)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_brewer(palette="Set1") +
  xlab("Provenance") + ylab("DBB (mm)")
OriginProvenancePlot_DBB

### Extra Trial 

###Height 

my_comparisons <- list( c("Inveraray", "New Forest"), c("Inveraray", "Rheidol"), c("New Forest", "Rheidol") )

CollectionSitePlot_Height <- ggplot(GlasshouseTrial_TheterophyllaExtra, aes(x=GlasshouseTrial_TheterophyllaExtra$Collection.site..UK., y=GlasshouseTrial_TheterophyllaExtra$Height..mm., fill=GlasshouseTrial_TheterophyllaExtra$Collection.site..UK.)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Collection Site (UK)") + ylab("Height (mm)") +
  #stat_compare_means(method= "t.test", comparisons = my_comparisons, label.y = c(45, 50, 55), hide.ns = TRUE) +
  stat_compare_means(method= "anova", label.y = 60)
CollectionSitePlot_Height

##DBB

CollectionSitePlot_DBB <- ggplot(GlasshouseTrial_TheterophyllaExtra, aes(x=GlasshouseTrial_TheterophyllaExtra$Collection.site..UK., y=GlasshouseTrial_TheterophyllaExtra$DBB..mm., fill=GlasshouseTrial_TheterophyllaExtra$Collection.site..UK.)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Collection Site (UK)") + ylab("DBB (mm)") +
 # stat_compare_means(method= "t.test", comparisons = my_comparisons, label.y = c(0.75, 0.8, 0.85), hide.ns = TRUE) +
  stat_compare_means(method= "anova", label.y = 0.9)
CollectionSitePlot_DBB
