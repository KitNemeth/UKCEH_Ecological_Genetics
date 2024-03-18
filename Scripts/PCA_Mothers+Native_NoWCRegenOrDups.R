# Install and load necessary libraries if not already installed
# You can install these libraries with install.packages("package_name") if needed.

# Install and load the 'plotly' library for interactive plotting
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}
library(plotly)

# Install and load the 'ggfortify' library for data fortification in ggplot2
if (!requireNamespace("ggfortify", quietly = TRUE)) {
  install.packages("ggfortify")
}
library(ggfortify)

# Install and load the 'data.table' library for data manipulation
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

# Install and load the 'svglite' library for generating SVG output
if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}
library(svglite)

# Set the working directory to your desired path
setwd("C:/Users/krinem/OneDrive - UKCEH")

# Read in TASSEL 5 PCA values in Plink format from a text file into a data frame
df <- read.table("PCA_Mothers+Native_NoWCRegenOrDups.txt", header = TRUE, sep = "\t")


# Read a Keyfile containing Taxa names and corresponding information of interest into a data frame
Keyfile <- read.table("newLEAFgenotypes_Keyfile_Mothers+Native.txt", header = TRUE, sep = "\t")

# Merge the two data frames based on specific columns
merged_data <- merge(df, Keyfile, by.x = "FID", by.y = "Taxa", row.names = FALSE)
merged_data$Source2 <- factor(merged_data$Source2, levels=c("Scotland", "Beanley Moss", "Sheringham Estate", "Felbrigg Estate", "Northern Scottish population", "Williams Cleugh"))


williams_cleugh_data <- subset(merged_data, Source2 == "Williams Cleugh")

# Create a ggplot2 plot
p <- ggplot(merged_data, aes(x = PC1, y = PC2, color = Source2)) +
  geom_point() +
  labs(x = "Principal Component 1 (1.0240%)",
       y = "Principal Component 2 (0.0996%)",
       title = "PCA Results") +
  scale_color_manual(values = c("#b69140",
                                 "#8d70c9",
                                 "#6aa74d",
                                 "#c8588c",
                                 "#49adad",
                                 "#cc5b43")) +
  stat_ellipse(data = williams_cleugh_data) +
  labs(color = "Source")  # Change the legend title to "Source"

# Convert the ggplot2 plot to an interactive Plotly plot
interactive_plot <- ggplotly(p, originalData = FALSE)  # Set originalData to FALSE

# Display the interactive plot
interactive_plot

# Define the output file path for the SVG file
out = "PCA_Scotland_Proportional_NoWilliamsCleughRegen.svg"

# Generate an SVG file from the ggplot2 plot
svg(out)

p

# Close the SVG device
dev.off()
