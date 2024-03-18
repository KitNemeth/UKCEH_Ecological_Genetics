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
setwd("C:/Users/krinem/OneDrive - UKCEH/PCA_Europe+Scotland_Proportional_NoWilliamsCleughRegen")

# Read in TASSEL 5 PCA values in Plink format from a text file into a data frame
df <- read.table("Europe+Scotland_Proportional_PCAValues_NoWilliamsCleughRegen.txt", header = TRUE, sep = "\t")

setwd("C:/Users/krinem/OneDrive - UKCEH")


# Read a Keyfile containing Taxa names and corresponding information of interest into a data frame
Keyfile <- read.table("Europe+Scotland_Proportional_Keyfile.txt", header = TRUE, sep = "\t")

# Merge the two data frames based on specific columns
merged_data <- merge(df, Keyfile, by.x = "FID", by.y = "Taxa", row.names = FALSE)

merged_data$Source <- factor(merged_data$Source, levels=c("Ireland", "Scotland", "Beanley Moss", "Williams Cleugh", "Sheringham Estate", "Felbrigg Estate", "Northern Scottish population", "Spain", "Italy", "Switzerland", "Austria", "Poland", "Sweden", "Finland", "Siberia", "Turkey"))

williams_cleugh_data <- subset(merged_data, Source == "Williams Cleugh")

# Create a ggplot2 plot
p <- ggplot(merged_data, aes(x = PC1, y = PC2, color = Source)) +
  geom_point() +
  scale_color_manual(values = c("#4aa98c",
                                "#ae58c7",
                                "#75b63b",
                                "#d9478a",
                                "#58bb6c",
                                "#636fc7",
                                "#b7ae49",
                                "#984f8a",
                                "#607c34",
                                "#cf8dd0",
                                "#df9445",
                                "#5ba5d7",
                                "#d14a3a",
                                "#9c6832",
                                "#a74557",
                                "#e38283")) +
  labs(x = "Principal Component 1 (3.42%)",
       y = "Principal Component 2 (2.44%)",
       title = "PCA Results") +
  stat_ellipse(data = williams_cleugh_data)

# Convert the ggplot2 plot to an interactive Plotly plot
interactive_plot <- ggplotly(p, originalData = FALSE)  # Set originalData to FALSE

# Display the interactive plot
interactive_plot

# Define the output file path for the SVG file
out = "PCA_Europe+Scotland_Proportional_NoWilliamsCleughRegen.svg"

# Generate an SVG file from the ggplot2 plot
svg(out)

p

# Close the SVG device
dev.off()
