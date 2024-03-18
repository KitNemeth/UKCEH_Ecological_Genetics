# Load the 'data.table' library for data manipulation
library(data.table)

# Read data from a file named "N1131_20795SNPs_renamed.hmp" into the 'GapGenotypes' data table
GapGenotypes <- fread("N1131_20795SNPs_renamed.hmp", header = TRUE, sep = "\t")

# Read data from a file named "newLEAFgenotypes.hmp.txt" into the 'newLEAFGenotypes' data table
newLEAFGenotypes <- fread("newLEAFgenotypes.hmp.txt", header = TRUE, sep = "\t")

# Merge the 'GapGenotypes' and 'newLEAFGenotypes' data tables by a common column "rs#"
merged_genotypes <- merge(GapGenotypes, newLEAFGenotypes, by = "rs#")

# Note: The next three lines are commented out (prefixed with #) and not active in the script.

# Uncommenting these lines would merge the tables with different suffixes for overlapping columns.
# merged_genotypes <- merge(GapGenotypes, newLEAFGenotypes, by = "rs#", suffixes = c("", ""))

# Remove specific columns (e.g., 'alleles.y', 'chrom.y', etc.) from the merged data table
# merged_genotypes <- merged_genotypes[, -c('alleles.y', 'chrom.y', 'pos.y', 'strand.y', 'assembly#.y', 'center.y', 'protLSID.y', 'assayLSID.y', 'panelLSID.y', 'Qccode.y')]

# Get the current column names of the merged data table
current_names <- names(merged_genotypes)

# Remove the ".x" suffix from the column names for those that have it
new_names <- sub(".x$", "", current_names)

# Assign the new column names to the data frame
names(merged_genotypes) <- new_names

# Sort the merged data table based on the "pos.x" column
merged_genotypes <- merged_genotypes[order(merged_genotypes$pos.x), ]

# Write the merged data table to a new file named "mergedgenotypes.hmp.txt" with tab-separated values, excluding row names and quotes
write.table(merged_genotypes, file = "mergedgenotypes.hmp.txt", sep = "\t", row.names = FALSE, quote = FALSE)