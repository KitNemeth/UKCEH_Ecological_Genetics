library(pheatmap)
library(data.table)

matrix <- as.matrix(fread("KinshipMatrix_Mothers+NativePine.txt"),rownames=1)

p <- pheatmap(matrix, display_numbers = F)


colnames(matrix) <- rownames(matrix)

order_list <- fread("Kinship_TaxaList_Order2.txt", header = FALSE, sep = "\t")[, 1]
order_list <- order_list$V1

# Reorder columns based on desired order
matrix_reordered <- matrix[, order(match(colnames(matrix), order_list))]
matrix_reordered2<- matrix_reordered[match(order_list, rownames(matrix_reordered)), ]

# Print the reordered matrix
print(matrix_reordered)

library(ASRgenomics)
p <- kinship.heatmap(
  K = matrix_reordered2,
  dendrogram = FALSE,
  clustering.method = c("hierarchical"),
  dist.method = c("euclidean"),
  row.label = TRUE,
  col.label = TRUE
)

matrix_reordered2 <-as.data.frame(matrix_reordered2)
 
rownames(matrix_reordered2) <- colnames(matrix_reordered2)

write.table(matrix_reordered2,"Kinship_Reordered_Mother+NativePine.txt",sep="\t",row.names=TRUE)
