library(ggplot2)
library(ggmap)
library(svglite)

UnitedKingdom <- c(left = -7.5, bottom = 50, right = 2, top = 59)
UnitedKingdomMap <- get_stamenmap(UnitedKingdom, zoom = 7, maptype = "terrain")
ggmap(UnitedKingdomMap)
  

p <- ggmap(UnitedKingdomMap) +
  geom_point(data = Keyfile,size=1.5, stat = "unique", shape=16, alpha=1, aes(x = Long, y = Lat, colour = combinedLatLong_color)) +  
  scale_colour_identity() 
p
ggsave(file="Map_LatLongColour.svg", plot=p, width = 6, height = 6)
