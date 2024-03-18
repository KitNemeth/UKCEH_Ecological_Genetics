library(plotly)
library(ggfortify)
library(data.table)
library(svglite)

df <- read.table("Europe+Scotland_Proportional_PCAValues.txt", header=T, sep="\t")

Keyfile <- read.table("Europe+Scotland_Proportional_Keyfile.txt", header=T, sep="\t")

merged_data <- merge(df, Keyfile, by.x = "FID", by.y = "Taxa", row.names = NULL)

p <- ggplot(merged_data, aes(x = PC1, y = PC2, color = Source)) +
  geom_point() +
  labs(x = "Principal Component 1",
       y = "Principal Component 2",
       title = "PCA Results")
interactive_plot <- ggplotly(p, originalData = FALSE)  # Set originalData to FALSE
interactive_plot
