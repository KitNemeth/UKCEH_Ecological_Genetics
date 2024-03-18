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

# Read in TASSEL 5 PCA values in Plink format from a text file into a data frame
df <- read.table("Europe+Scotland_PCAValues.txt", header = TRUE, sep = "\t")

# Read a Keyfile containing Taxa names and corresponding information of interest into a data frame
Keyfile <- read.table("Europe+Scotland_Keyfile.txt", header = TRUE, sep = "\t")

# Merge the two data frames based on specific columns
merged_data <- merge(df, Keyfile, by.x = "FID", by.y = "Taxa", row.names = FALSE)

# Create a ggplot2 plot
p <- ggplot(merged_data, aes(x = PC1, y = PC2, color = Source)) +
  geom_point() +
  labs(x = "Principal Component 1",
       y = "Principal Component 2",
       title = "PCA Results")

# Convert the ggplot2 plot to an interactive Plotly plot
interactive_plot <- ggplotly(p, originalData = FALSE)  # Set originalData to FALSE

# Display the interactive plot
interactive_plot

# Define the output file path for the SVG file
out = "PCA_Europe+Scotland.svg"

# Generate an SVG file from the ggplot2 plot
svg(out)

p

# Close the SVG device
dev.off()