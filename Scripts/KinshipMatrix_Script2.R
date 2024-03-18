data <- read.table("PineTrial_KinshipMatrix.txt", skip = 3)

kinship_matrix <- as.matrix(data[1:nrow(data), -c(1)])


# Create the heatmap
heatmap(
  as.matrix(kinship_matrix),  # Convert to matrix
  col = colorRampPalette(c("white", "blue"))(256),  # Define the color palette
  scale = "none",  # Use "none" if you don't want to scale the values
  main = "Kinship Matrix Heatmap",  # Set the title of the heatmap
  xlab = "Sample IDs",  # Label for x-axis
  ylab = "Sample IDs"   # Label for y-axis
)
