GlasshouseTrial_Psylvestris <- read.delim("C:/Users/krinem/OneDrive - UKCEH/GlasshouseTrial_Psylvestris.txt")

library(ggpubr)
library(ggplot2)

GlasshouseTrial_Psylvestris$Soil <- factor(GlasshouseTrial_Psylvestris$Soil , levels=c("L", "M", "H"))

my_comparisons <- list( c("L", "M"), c("M", "H"), c("L", "H") )

GlasshouseTrial_Psylvestris$Height.adjusted..mm. <- as.numeric(GlasshouseTrial_Psylvestris$Height.adjusted..mm.)


HeightPlot <- ggplot(GlasshouseTrial_Psylvestris, aes(x=Soil, y=Height.adjusted..mm., fill=Soil)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Soil Fertility") + ylab("Height (mm)") +
  stat_compare_means(method= "t.test", comparisons = my_comparisons, label.y = c(80, 85, 90)) +
  stat_compare_means(method= "anova", label.y = 100)
HeightPlot

GlasshouseTrial_Psylvestris$DBB..mm..01.02.24 <- as.numeric(GlasshouseTrial_Psylvestris$DBB..mm..01.02.24)

DiameterPlot <- ggplot(GlasshouseTrial_Psylvestris, aes(x=Soil, y=DBB..mm..01.02.24, fill=Soil)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  xlab("Soil Fertility") + ylab("DBB (mm)") +
  stat_compare_means(method= "t.test", comparisons = my_comparisons, label.y = c(1.3, 1.35, 1.4)) +
  stat_compare_means(method= "anova", label.y = 1.5)
DiameterPlot

