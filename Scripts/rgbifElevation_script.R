library(rgbif)
library(openxlsx)

Elev <- elevation(latitude = Keyfile$Lat,
          longitude = Keyfile$Long,
          username = "krinem")

Keyfile <- cbind(Keyfile, Elev = Elev$elevation_geonames)

write.xlsx(Keyfile, 'newLEAFgenotypes_Keyfile.xlsx')
