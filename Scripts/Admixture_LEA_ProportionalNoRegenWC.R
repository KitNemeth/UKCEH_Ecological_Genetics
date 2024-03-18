# Install the BiocManager package if it's not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the LEA package from Bioconductor
BiocManager::install("LEA")

# Load the LEA library
library(LEA)

# Set the working directory
setwd("C:/Users/krinem/OneDrive - UKCEH")

# Convert VCF file to LFMM format and save it
vcf2lfmm("Europe+Scotland_FilteredGenotypes_Proportional_NoWCRegen.vcf", output.file = "Europe+Scotland_FilteredGenotypes_NoRegenWC", force = TRUE)

# Perform Principal Component Analysis (PCA) on LFMM data
pc = pca("Europe+Scotland_FilteredGenotypes_NoRegenWC.lfmm", scale = TRUE)

# Calculate Tracy-Widom statistics
tw = tracy.widom(pc)

# Plot the percentage of variance explained by each component
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)

# Create an SNMF project
project = snmf("Europe+Scotland_FilteredGenotypes_NoRegenWC.geno",
               K = 1:12,
               entropy = TRUE,
               repetitions = 10,
               project = "new")

# Save a plot of the cross-entropy criterion
out = "cross-entropy criterion.svg"
svg(out)

# Plot the cross-entropy criterion for all runs in the SNMF project
plot(project, col = "blue", pch = 19, cex = 1.2)

# Turn off the SVG device
dev.off()

# Select the best run for K = 7 clusters
best = which.min(cross.entropy(project, K = 7))

# Read population data from a file
Populations_data <- read.table("geno_poplist_NoRegenWC.txt", header = TRUE, sep = "\t")

# Read the desired order from a file
desired_order <- read.csv("desired order_NoRegenWC.txt", header = FALSE)[, 1]

# Calculate the Q-matrix for the best run
qmatrix <- LEA::Q(project, K = 7, run = best)
qmatrix <- as.data.frame(qmatrix)

# Combine population data and Q-matrix, and sort by desired order
tbl <- cbind(Populations_data, qmatrix)
tbl <- tbl[order(match(tbl$taxa, desired_order)), ]

# Define colors for the barplot
col = c("#48b1a7", "#b05cc6", "#6ca74d", "#6f7ccb", "#b68f40", "#c65c8a", "#cc5a43")

# Find breakpoints between different populations
breaks <- 0
for (i in 2:nrow(tbl)) {
  if (tbl$pop[i] != tbl$pop[i - 1]) {
    breaks <- c(breaks, i - 1)
  }
}
breaks <- c(breaks, nrow(tbl))

# Calculate spaces and labels for the barplot
spaces <- rep(0, nrow(tbl))
spaces[breaks[2:(length(breaks) - 1)] + 1] <- 2
labels <- rep("", nrow(tbl))
labels[breaks[1:(length(breaks) - 1)] + ((breaks[2:length(breaks)] - breaks[1:(length(breaks) - 1)]) / 2)] <- gsub("_", " ", tbl$pop[breaks])

# Save the barplot as an SVG
out = "Admixture_Europe+Scotland_FilteredGenotypes_NoRegenWC.svg"
svg(out, width = 20, height = 10)

# Set plot parameters
par(oma = c(4.5, 1, 1, 1), mar = c(7, 3, 0, 7), mgp = c(2, 1, 0), xpd = TRUE)

# Create the barplot
barplot(t(as.matrix(tbl[, 3:9])), col = col, xaxt = "n", border = NA, xlim = c(0, nrow(tbl)), width = 1, space = spaces, names.arg = rep("", nrow(tbl)), ylab = "Admixture coefficients")

# Add x-axis labels with line breaks
axis(1, line = -0.2, at = breaks[1:(length(breaks) - 1)] + ((0:(length(breaks) - 2)) * 2) + ((breaks[2:length(breaks)] - breaks[1:(length(breaks) - 1)]) / 2), lwd = 0, lwd.ticks = 1, labels = gsub("_", " ", tbl$pop[breaks]), las = 2)

# Turn off the SVG device
dev.off()