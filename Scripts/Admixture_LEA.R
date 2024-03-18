install.packages(c("fields","RColorBrewer","mapplots"))


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

browseVignettes("LEA")

library(LEA)

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/PopulationsSutilities.R")

library(pegas)

install.packages("VariantAnnotation")
library(VariantAnnotation)

setwd("C:/Users/krinem/OneDrive - UKCEH")

vcf2lfmm("Europe+Scotland_FilteredGenotypes_Proportional.vcf", output.file = "Europe+Scotland_FilteredGenotypes_Proportional", force = TRUE)

pc = pca("Europe+Scotland_FilteredGenotypes_Proportional.lfmm", scale = TRUE)
tw = tracy.widom(pc)

#plot the percentage of variance explained by each component
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)

project = NULL
project = snmf("Europe+Scotland_FilteredGenotypes_Proportional.geno",
               K = 1:12,
               entropy = TRUE,
               repetitions = 10,
               project = "new")

out="cross-entropy criterion.svg"
svg(out)
# plot cross-entropy criterion for all runs in the snmf project
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()

# select the best run for K = 4 clusters
best = which.min(cross.entropy(project, K = 7))
my.colors <- c("#cd5760",
               "#89aa3c",
               "#8a67d2",
               "#5fa271",
               "#c55c9c",
               "#c2753a",
               "#718cc7")

barchart(project, K = 7, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") 
axis(1, at = 1:length(my_list),
     labels = my_list, las=1,
     cex.axis = .4)

Populations_data <-  read.table("geno_poplist.txt", header = TRUE, sep = "\t")
desired_order <- read.csv("desired order.txt", header = FALSE)[, 1]

qmatrix = LEA::Q(project, K = 7, run = best)
qmatrix <-as.data.frame(qmatrix)
tbl <- cbind(Populations_data,qmatrix)
tbl <- tbl[order(match(tbl$taxa, desired_order)), ]
tbl <- as.data.frame(tbl)


pops <- unique(unlist(lapply(tbl$pop, function(x) {
  y <- c()
  y <- c(y, unlist(strsplit(x, "_")[[1]][1]))
})))       

sample_sites <- rep(NA,nrow(tbl))
  for (i in 1:nrow(tbl)){
  sample_sites[i] <- strsplit(rownames(tbl),"_")[[i]][1]}
N <- unlist(lapply(pops,function(x){length(which(sample_sites==x))}))
names(N) <- pops
N


# Set names for N based on pops
names(N) <- pops
N

#par(mar=c(4,4,0.5,0.5))

tbl$pop <- factor(tbl$pop)


col= c("#48b1a7",
       "#b05cc6",
       "#6ca74d",
       "#6f7ccb",
       "#b68f40",
       "#c65c8a",
       "#cc5a43")

k <- 7

breaks <- 0
for (i in 2:nrow(tbl)) {
  if (tbl$pop[i]!=tbl$pop[i-1]) {
    breaks <- c(breaks,i-1)
  }
}
breaks <- c(breaks,nrow(tbl))
#breaks <- 0
#spaces <- 0

spaces <- rep(0,nrow(tbl))
spaces[breaks[2:(length(breaks)-1)]+1] <- 2
labels <- rep("",nrow(tbl))
labels[breaks[1:(length(breaks)-1)]+((breaks[2:length(breaks)]-breaks[1:(length(breaks)-1)])/2)] <- gsub("_"," ",tbl$pop[breaks])


out="Admixture_Europe+Scotland_FilteredGenotypes_Proportional.svg"
svg(out,width=20,height=10)

par(oma = c(4.5, 1, 1, 1),   # outer margins for titles and axis labels
    mar = c(7, 3, 0, 7),   # margins for the bottom, left, top, and right sides
    mgp = c(2, 1, 0),      # axis label at 2 rows distance, tick labels at 1 row
    xpd = TRUE)            # allow content to protrude into outer margin (and beyond)
barplot(t(as.matrix(tbl[, 3:9])), col=col, xaxt="n",border=NA,xlim=c(0,nrow(tbl)),width = 1,space = spaces,names.arg=rep("",nrow(tbl)), ylab="Admixture coefficients")        
axis(1,line=-0.2,at = breaks[1:(length(breaks)-1)]+((0:(length(breaks)-2))*2)+((breaks[2:length(breaks)]-breaks[1:(length(breaks)-1)])/2),lwd=0,lwd.ticks=1,labels = gsub("_"," ",tbl$pop[breaks]), las=2)

dev.off()



