install.packages("vcfR")
library(vcfR) #using the package vcfR, which allows you to view and analyze vcf files
vcf <- read.vcfR("Europe+Scotland_FilteredGenotypes_Proportional.vcf")
PS_genind <- vcfR2genind(vcf)
PS_genind
PS_genlight <- vcfR2genlight(vcf) # more recently, adegenet package uses genlight objects for highthroughput analysis such as SNP's.

install.packages("adegenet")
library(adegenet)

# Run the admixture analysis
admixture_result <- dapc(PS_genlight, n.pca = 5, n.da = 3)  # You can adjust the parameters

# Visualize the admixture results
plot(admixture_result)

library(writexl)
taxa_list <- as.data.frame(PS_genlight$ind.names)
write_xlsx(taxa_list, path = "taxa_list2.xlsx")

pop_data <- read.table("pop.txt", header = FALSE, sep = "\t", col.names = "Population")
pop_vector <- pop_data$Population

pop(PS_genlight) <- pop_vector

grp <- find.clusters(PS_genind, max.n.clust=40)

dapc1 <- dapc(PS_genind, grp$grp)

scatter(dapc1)

myCol <- c("darkblue","purple","green","orange","red","blue")


scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=myCol, solid=.4, cex=3, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, txt.leg=paste("Cluster",1:6))
par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=myCol)

compoplot(dapc1, posi="bottomright",
          txt.leg=paste("Cluster", 1:6), lab="",
          ncol=1, xlab="individuals", col=funky(6))
