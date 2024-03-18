library(plotly)
library(ggfortify)
library(data.table)
library(svglite)

df <- read.table("PCA_Values_newLEAF_MothersNative.txt", header=T, sep="\t")

df[2] <- NULL 
df[1] <- NULL 

pca_res <- prcomp(df, scale. = TRUE)

Keyfile <- fread("newLEAFgenotypes_Keyfile.txt", header=T, sep="\t")
rownames(Keyfile) <- Keyfile[,1]


r <- ((Keyfile[,Lat]-min(Keyfile[,Lat]))/(max(Keyfile[,Lat])-min(Keyfile[,Lat])))*255
g <- ((Keyfile[,Long]-min(Keyfile[,Long]))/(max(Keyfile[,Long])-min(Keyfile[,Long])))*255
b <- rep(0,nrow(Keyfile))
print(rgb(r,g,b,maxColorValue = 255))
return(rgb(r,g,b,maxColorValue = 255))
Keyfile$combinedLatLong_color <- (rgb(r,g,b,maxColorValue = 255))
View(Keyfile)

p <- autoplot(pca_res, data = Keyfile, colour = 'combinedLatLong_color') + 
  labs(x = "PC1 (1.35%)", y = "PC2 (1.24%)") +  scale_colour_identity()

out="PCA_LatLongColour.svg"
svg(out)

PCA<-ggplotly(p) 
PCA

dev.off()

ggsave(file="PCA_LatLongColour.svg", plot=p)

p <- autoplot(pca_res, data = Keyfile, colour = 'Source') + 
  labs(x = "PC1 (1.35%)", y = "PC2 (1.24%)")
scale_colour_identity()


merged_data <- merge(df, Keyfile, by.x = "FID", by.y = "Taxa", row.names = NULL)
rownames(merged_data) <- df[, 1]
merged_data[1] <- NULL 
p <- ggplot(merged_data, aes(x = PC1, y = PC2, color = merged_data$Source)) +
  geom_point() +
  labs(x = "Principal Component 1",
       y = "Principal Component 2",
       title = "PCA Results")

df <- read.table("PCA_Values_newLEAF_All.txt", header=T, sep="\t")

