source("C:/Users/krinem/OneDrive - UKCEH/Scripts/POPSutilities.R")

install.packages(c("fields","RColorBrewer","mapplots"))
library(fields)
library(RColorBrewer)
library(mapplots)
library(maps)

library(devtools)
install.packages("devtools")
devtools::install_github("bcm-uga/TESS3_encho_sen")

install.packages("remotes")
remotes::install_github("cayek/TESS3_encho_sen")

# Load the tess3r library
library(tess3r)

# Install the LEA package from Bioconductor
BiocManager::install("LEA")

# Load the LEA library
library(LEA)
setwd("C:/Users/krinem/OneDrive - UKCEH/Admixture map LEA")

project = load.snmfProject("Europe+Scotland_FilteredGenotypes_Proportional_NoWCRegen.snmfProject")
best = which.min(cross.entropy(project, K = 6))

qmatrix <- LEA::Q(project, K = 6, run = best)

asc.raster="http://membres-timc.imag.fr/Olivier.Francois/RasterMaps/Europe.asc"
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster,cell_value_min=0)

coord = read.table("Test.coords.txt") 
coord$V1 <- as.numeric(coord$V1)
coord$V2 <- as.numeric(coord$V2)


maps(matrix = qmatrix, coord, grid, method = "max", main = "Ancestry coefficients", xlab = "Longitude (°E)", ylab = "Latitude (°N)", cex = .5)
map(add = T, interior = F)



show.key = function(cluster=1,colorGradientsList=lColorGradients){
  ncolors=length(colorGradientsList[[cluster]])
  barplot(matrix(rep(1/10,10)),col=colorGradientsList[[cluster]][(ncolors-9):ncolors],main=paste("Cluster",cluster))}

# Increase the size of the plotting device
svg("AncestryMap.svg", width = 10, height = 8)

# Adjust the layout to ensure proper margins
layout(
  matrix(c(rep(1, 12), 2, 3, 4, 5, 6, 7), 6, 3, byrow = FALSE),
  widths = c(3, 1),
  respect = TRUE
)

# Set plotting parameters
par(ps = 22, mar = c(5, 5, 2, 2))  # Adjust mar to set margins

# Your plotting code
maps(matrix = qmatrix, coord, grid, method = "max", main = "Ancestry coefficients", xlab = "Longitude (°E)", ylab = "Latitude (°N)", cex = 0.5)
map(add = TRUE, interior = FALSE)

par(ps = 16)

# Loop for legends
for (k in 1:6) {
  show.key(k)
}

dev.off()
# Reset the layout
layout(1)

# Close the plotting device
dev.off()
# Reset the layout
layout(1)

asc.raster <- tempfile()
download.file("http://membres-timc.imag.fr/Olivier.Francois/RasterMaps/Europe.asc", asc.raster)
plot(qmatrix, coord, method = "map.max", cex = .4,  
     interpol = FieldsKrigModel(10), 
     raster.filename = asc.raster,
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude", 
     col.palette = my.palette)
