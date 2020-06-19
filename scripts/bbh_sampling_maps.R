library(maps)
library(mapdata)
library(rgdal)

# Read in metadata for landlocked and baseline BBH
meta_land <- read.table("~/Documents/UCSC_postdoc/blueback_herring/data/landlocked_metadata.txt")
# meta_base <- read.delim("~/Documents/UCSC_postdoc/blueback_herring/data/baseline_metadata.txt", header = TRUE)
base_latlon <- read.delim("~/Documents/UCSC_postdoc/blueback_herring/maps/baseline_latlon.txt", header = TRUE)

# Read in shapefiles with river and lake data
riversData <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp") # load the shapefile
riversData_fine <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/ne_10m_rivers_north_america/ne_10m_rivers_north_america.shp")
lakesData <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/ne_10m_lakes/ne_10m_lakes.shp")
lakesData_fine <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/ne_10m_lakes_north_america/ne_10m_lakes_north_america.shp")

rs <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/Lakes%2C_Ponds%2C_Reservoirs%2C_and_Swamps_Georgia-shp/Lakes%2C_Ponds%2C_Reservoirs%2C_and_Swamps_Georgia.shp")



# Find unique lat/long combos
coors <- meta_land[!duplicated(meta_land[8:10]),]
coors <- coors[-c(5,11:12),] # if lakes have multiple lat/lon, keep most common: Nottley (34.93285, -84.03903); Tugalo (34.75695, -83.31985)
# coors_base <- meta_base[!duplicated(meta_base[5:7]),]

#### Large map ####
png(file="~/Documents/UCSC_postdoc/blueback_herring/maps/anadromous.png", width=9, height=7.5, res=300, units="in")

par(
  mar=c(5, 7, 4, 2), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

map("worldHires", c("us", "canada"), xlim=c(-86,-60), ylim=c(28,48), col="gray92", fill=TRUE) #plots the region of the USA that I want
map("state", xlim=c(-86,-60), ylim=c(28,48), add = TRUE, boundary=FALSE, col = 'gray70') # plots US state boundaries
plot(riversData, col='skyblue2', add=T) # plot big rivers
title(xlab = "Longitude (°)", ylab = "Latitude (°)")

axis(1, at=seq(-85,-60, by=5), labels=seq(-85,-60, by= 5))
axis(2, at=seq(28,48, by = 5), labels=seq(28,48, by= 5), las = TRUE)
box()

# Plot sampling locations
rect(-84.5,34,-82.5,35.1, border = 'tomato', lwd = 1.5) # box around inset
lines(c(-70.91866, -70), c(43.13086,43.3)) # lines first so the points are plotted on top
lines(c(-70.94437, -70), c(42.98110,42.9))
lines(c(-70.92917, -70), c(42.75007,42.4))
lines(c(-69.5, -70.12213), c(41.7, 41.68207))
lines(c(-70.62445,-70.6),c(41.73676,40.8))
points(base_latlon$Longitude, base_latlon$Latitude, col = 'black', pch = 21, bg = 'tomato')
points(coors$Longitude,coors$Latitude, col = 'black', pch = 21, bg = 'tomato')

# Add labels
text(-65,40, expression(italic("Atlantic Ocean")), cex = 0.8)
text(-61.8, 46.8,"MAR", cex = 0.7)
text(-65.5, 46.3,"PET", cex = 0.7)
text(-68.25, 44.9,"EMA", cex = 0.7)
text(-70.5, 44.5,"LOC", cex = 0.7)
text(-69.3, 43.4,"OYS", cex = 0.7)
text(-69.3, 42.9,"EXE", cex = 0.7)
text(-69.3, 42.4,"PAR", cex = 0.7)
text(-72.1, 42.4,"MYS", cex = 0.7)
text(-68.8, 41.7,"HER", cex = 0.7)
text(-70.6, 40.5,"MON", cex = 0.7)
text(-72.2, 41.7,"GIL", cex = 0.7)
text(-73.2, 40.2,"MET", cex = 0.7)
text(-75.7, 40.15,"DEL", cex = 0.7)
text(-77, 39.7,"SUS", cex = 0.7)
text(-78, 39,"POT", cex = 0.7)
text(-78.3, 34.6,"CPF", cex = 0.7)
text(-80, 33.7,"SAN", cex = 0.7)
text(-81.5, 32.8,"SAV", cex = 0.7)
text(-82.1, 30.1,"STJ", cex = 0.7)

dev.off()

#### Inset map ####
map("worldHires", c("us"), xlim=c(-84.5,-82.5), ylim=c(34,35.1), col="gray90", fill=TRUE) #plots the region of the USA that I want
map("stdate", xlim=c(-84.5,-82.5), ylim=c(34,35.1), add = TRUE, boundary=FALSE, col = 'black') # plots US state boundaries
title(xlab = "Longitude (°)", ylab = "Latitude (°)")

axis(1, at=seq(-85,-82, by=1), labels=seq(-85,-82, by= 1))
axis(2, at=seq(33,36, by = 1), labels=seq(33,36, by= 1), las = TRUE)
box()

# Plot rivers and lakes
plot(riversData, col='blue', add=T) # plot big rivers
plot(riversData_fine, col='blue', add=T) # plot smaller rivers
plot(lakesData, col='blue', add=T)
plot(lakesData_fine, col='blue', add=T)

plot(rs, col = 'blue', add=T)

# Plot sampling locations
points(coors$Longitude, coors$Latitude, pch= 19, col= 'tomato')

