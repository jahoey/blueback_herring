library(maps)
library(mapdata)
library(rgdal)

#### Read in data ####
# Read in metadata for landlocked and baseline BBH
meta_land <- read.table("~/Documents/UCSC_postdoc/blueback_herring/data/landlocked_metadata.txt")
# meta_base <- read.delim("~/Documents/UCSC_postdoc/blueback_herring/data/baseline_metadata.txt", header = TRUE)
anadromous_latlon <- read.delim("~/Documents/UCSC_postdoc/blueback_herring/maps/anadromous_latlon.txt", header = TRUE) # Many anadromous populations from Reid et al 2018, but not all

# Read in metadata for landlocked and baseline ALE
ale_base_latlon <- read.delim("~/Documents/UCSC_postdoc/alewife/data/ALE_anadromous_latlon.txt", header = TRUE, fileEncoding="latin1")
ale_land_latlon <- read.delim("~/Documents/UCSC_postdoc/alewife/data/ALE_landlocked_latlon.txt", header = TRUE)

# Which sites are shared by BBH and ALE?
both <- intersect(ale_base_latlon$Location, anadromous_latlon$Code)

# Make different data frames for BBH only, ALE only and both
anadromous_latlon_only <- anadromous_latlon[!anadromous_latlon$Code %in% both,] # BBH only
ale_base_latlon_only <- ale_base_latlon[!ale_base_latlon$Location %in% both,] # ALE only
both_only <- ale_base_latlon[ale_base_latlon$Location %in% both,] # both BBH & ALE

# Read in shapefiles with river and lake data
us_rivers <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp")
us_rivers_fine <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/ne_10m_rivers_north_america/ne_10m_rivers_north_america.shp")
H_0602Lakes <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/H_0602_lakes/H_0602_lakes.shp")
H_0306Lakes <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/H_0306_lakes/H_0306_lakes.shp")
H_0313Lakes <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/H_0313_lakes/H_0313_lakes.shp")
H_0315Lakes <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/H_0315_lakes/H_0315.shp")
H_0307Lakes <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/H_0307_lakes/H_0307.shp")
H_0305Lakes <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/H_0305_lakes/H_0305.shp")
H_0601Lakes <- readOGR("~/Documents/UCSC_postdoc/blueback_herring/maps/H_0601_lakes/H_0601.shp")

# Find unique lat/long combos for BBH
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

map("worldHires", c("us", "canada"), xlim=c(-90,-60), ylim=c(28,48), col="gray92", fill=TRUE) #plots the region of the USA that I want
map("state", xlim=c(-90,-60), ylim=c(28,48), add = TRUE, boundary=FALSE, col = 'gray70') # plots US state boundaries
plot(us_rivers, col='skyblue2', add=T) # plot big rivers
# plot(us_rivers_fine, col='skyblue2', add=T)
title(xlab = "Longitude (°)", ylab = "Latitude (°)")

axis(1, at=seq(-90,-60, by=5), labels=seq(-90,-60, by= 5))
axis(2, at=seq(28,48, by = 5), labels=seq(28,48, by= 5), las = TRUE)
box()

# Plot BBH sampling locations
rect(-84.5,34,-82.5,35.1, border = 'black', lwd = 1.5) # box around inset
rect(-73,40.7,-71.9,41.6, border = 'black', lwd = 1.5) # box around inset
lines(c(-70.91866, -70), c(43.13086,43.3)) # lines first so the points are plotted on top
lines(c(-70.94437, -70), c(42.98110,42.9))
lines(c(-70.92917, -70), c(42.75007,42.4))
lines(c(-69.5, -70.12213), c(41.7, 41.68207))
# lines(c(-70.62445,-70.6),c(41.73676,40.8))
lines(c(-77.8, -76.8985), c(37.1,37.3727 ))
lines(c(-78.1, -76.7903), c(37.55,37.5161 ))
lines(c(-75.2, -76.13680), c(35.9, 35.69768))
lines(c(-77.9, -76.67790), c(38.5, 38.58754))
lines(c(-71.44630, -71), c(41.51840, 40.75))
lines(c(-70.62450, -70.2), c(41.73680, 40.95))
lines(c(-71.5, -71.13810), c(43, 42.41530))
lines(c(-71.6, -70.66420), c(42.3, 41.95425)) # TOW1
lines(c(-71.6, -71.77), c(42.3, 42.58)) # TOW2
lines(c(-71.08463, -72), c(41.89072, 42.2)) # NEM
lines(c(-75.7, -76.27955), c(43, 42.85962)) # OSL
lines(c(-77.28993, -78.1), c(42.79359, 42.4)) # CAL
lines(c(-77.6, -76.92872), c(41.95,42.69501)) #SNL
lines(c(-76, -76.73636), c(41.95,42.72954)) #CYL
lines(c(-70.3,-69.27661), c(44.35,44.23834)) # STG
lines(c(-67.2, -68.74302), c(42.4, 44.56971)) # ORL
lines(c(-67.1,-68.42883), c(43.1,44.54369)) # UNI
lines(c(-87,-87.85095), c(42.75,43.38235)) # LMI (left)
lines(c(-87,-86.49678), c(42.75,43.94479)) # LMI (right)
# Plot BBH lat/lon
points(anadromous_latlon_only$Longitude, anadromous_latlon_only$Latitude, col = 'black', pch = 21, bg = 'cadetblue3') # BBH anadromous only
points(coors$Longitude,coors$Latitude, col = 'black', pch = 21, bg = 'cadetblue3') # BBH landlocked
# Plot ALE lat/lon
points(ale_base_latlon_only$Longitude, ale_base_latlon_only$Latitude, col = 'black', pch = 21, bg = 'tomato') # ALE anadromous only
points(ale_land_latlon$Longitude, ale_land_latlon$Latitude, col = 'black', pch = 21, bg = 'tomato') # ALE landlocked
# Both
points(both_only$Longitude, both_only$Latitude, col = 'black', pch = 21, bg = 'slateblue2')

# Add labels
text(-65,40, expression(italic("Atlantic Ocean")), cex = 0.8)
text(-66.0, 47.35, "MIR", cex = 0.7)
text(-66.85, 46.4,"SJR", cex = 0.7)
text(-62.35, 45.4,"WAU", cex = 0.7)
text(-62.85, 44.45,"SUL", cex = 0.7)
text(-66.7, 45.4,"DEN", cex = 0.7)
text(-78.5, 38.5,"PAT", cex = 0.7)
text(-74.55, 35.9,"ALL", cex = 0.7)
text(-79.3, 36.6,"KER", cex = 0.7)
text(-68.6, 45.73, "EGL", cex = 0.7)
text(-74.85, 41.2, "HUD", cex =0.7)
text(-74.2, 44.54, "LCH", cex = 0.7)
text(-72.5, 42.7, "TOW", cex = 0.7)
text(-72.7, 42.2, "NEM", cex = 0.7)
text(-89.3, 47.45, "LSU", cex = 0.7)
text(-78.6, 43.6, "LON", cex = 0.7)
text(-75, 43, "OSL", cex = 0.7)
text(-78.6, 42.4, "CAL", cex = 0.7)
text(-77.6, 41.7, "SNL", cex = 0.7)
text(-76, 41.7, "CYL", cex = 0.7)
text(-70.95, 43.9, "AND", cex = 0.7)
text(-71, 44.4, "STG", cex = 0.7)
text(-66.5, 42.4, "ORL", cex = 0.7)
text(-66.5, 43, "UNI", cex = 0.7)
text(-87, 42.5, "LMI", cex = 0.7)

text(-61.8, 46.8,"MAR", cex = 0.7)
text(-65.5, 46.3,"PET", cex = 0.7)
text(-66.9, 44.4,"EMA", cex = 0.7)
# text(-70.5, 44.5,"LOC", cex = 0.7) # BBH only map
text(-69.9, 45,"LOC", cex = 0.7) # BBY and ALE map
text(-69.3, 43.4,"OYS", cex = 0.7)
text(-69.3, 42.9,"EXE", cex = 0.7)
text(-69.3, 42.4,"PAR", cex = 0.7)
text(-72.1, 43.2,"MYS", cex = 0.7)
text(-68.8, 41.7,"HER", cex = 0.7)
text(-70.6, 40.5,"GIL", cex = 0.7)
text(-69.4, 40.85,"MON", cex = 0.7)
text(-73.25, 39.85,"MET", cex = 0.7)
text(-75.7, 40.15,"DEL", cex = 0.7)
text(-77, 39.7,"SUS", cex = 0.7)
text(-78, 39,"POT", cex = 0.7)
text(-77.7, 38,"RAP", cex = 0.7)
text(-78.8, 37.55,"YOR", cex = 0.7)
text(-78.5, 37.1,"JAM", cex = 0.7)
text(-76.75, 36.55,"CHO", cex = 0.7)
text(-77.75, 36,"ROA", cex = 0.7)
text(-78.2, 35.3,"NEU", cex = 0.7)
text(-78.3, 34.6,"CF", cex = 0.7)
text(-80, 33.7,"SAN", cex = 0.7)
text(-81.5, 32.75,"SAV", cex = 0.7)
text(-82.4, 31.5,"ALT", cex = 0.7)
text(-82.1, 30.1,"STR", cex = 0.7)

# legend
legend("bottomright", c('Alewife', 'Both alewife & blueback herring', 'Blueback herring'), pch = c(21, 21, 21), pt.bg = c('tomato', 'slateblue2', 'cadetblue3'), cex = 0.8)

dev.off()

#### Inset map for southern reservoirs ####
png(file="~/Documents/UCSC_postdoc/blueback_herring/maps/landlocked.png", width=7, height=5, res=300, units="in")

par(
  mar=c(5, 7, 4, 2), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)
map("worldHires", c("us"), xlim=c(-84.5,-82.5), ylim=c(34,35.1), col="gray92", fill=TRUE) #plots the region of the USA that I want
map("state", xlim=c(-84.5,-82.5), ylim=c(34,35.1), add = TRUE, boundary=FALSE, col = 'gray70') # plots US state boundaries
title(xlab = "Longitude (°)", ylab = "Latitude (°)")

axis(1, at=seq(-85,-82, by=1), labels=seq(-85,-82, by= 1))
axis(2, at=seq(33,36, by = 1), labels=seq(33,36, by= 1), las = TRUE)
box(col = 'tomato')

# Plot rivers and lakes
plot(H_0602Lakes[which(H_0602Lakes$SHAPE_Area > 0.00001),], col = 'skyblue2', add = T, border = 'skyblue2')
plot(H_0306Lakes[which(H_0306Lakes$SHAPE_Area > 0.00001),], col = 'skyblue2', add = T, border = 'skyblue2')
plot(H_0313Lakes[which(H_0313Lakes$SHAPE_Area > 0.00001),], col = 'skyblue2', add = T, border = 'skyblue2')
plot(H_0315Lakes[which(H_0315Lakes$SHAPE_Area > 0.00001),], col = 'skyblue2', add = T, border = 'skyblue2')
plot(H_0307Lakes[which(H_0307Lakes$SHAPE_Area > 0.00001),], col = 'skyblue2', add = T, border = 'skyblue2')
plot(H_0305Lakes[which(H_0305Lakes$SHAPE_Area > 0.00001),], col = 'skyblue2', add = T, border = 'skyblue2')
plot(H_0601Lakes[which(H_0601Lakes$SHAPE_Area > 0.00001),], col = 'skyblue2', add = T, border = 'skyblue2')

# Plot sampling locations
points(coors$Longitude, coors$Latitude, pch= 21, bg= 'tomato', col = 'black')

# Add sampling labels
text(-83.925, 34.93,"LNO")
text(-83.638, 35,"LCT")
text(-83.68, 34.8,"LBU")
text(-83.445, 34.84,"LSE")
text(-83.49680, 34.7,"LRA")
text(-83.205, 34.761,"LTU")
text(-83.24, 34.68,"LYO")
text(-83.14807, 34.47,"LHA")
text(-83.8, 34.26,"LLA")

dev.off()

#### Inset map for Long Island Sound ####
png(file="~/Documents/UCSC_postdoc/blueback_herring/maps/li_sound.png", width=7, height=5, res=300, units="in")

par(
  mar=c(5, 7, 4, 2), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)
map("worldHires", c("us"), xlim=c(-72.9,-71.8), ylim=c(40.6,41.7), col="gray92", fill=TRUE) #plots the region of the USA that I want
map("state", xlim=c(-72.9,-71.8), ylim=c(40.6,41.7), add = TRUE, boundary=FALSE, col = 'gray70') # plots US state boundaries
title(xlab = "Longitude (°)", ylab = "Latitude (°)")

axis(1, at=seq(-72.9,-71.8, by=1), labels=seq(-72.9,-71.8, by= 1))
axis(2, at=seq(40.6,41.7, by = 1), labels=seq(40.6,41.7, by= 1), las = TRUE)
box(col = 'black')

# Plot rivers and lakes

# Plot sampling locations
points(ale_base_latlon_only$Longitude, ale_base_latlon_only$Latitude, col = 'black', pch = 21, bg = 'tomato') # ALE anadromous only
points(ale_land_latlon$Longitude, ale_land_latlon$Latitude, col = 'black', pch = 21, bg = 'tomato') # ALE landlocked

# Add sampling labels
text(-72.7, 41.45,"QUO")
text(-72.78, 40.9,"PEC")
text(-72.4, 41.41,"ROG")
text(-72.15, 41.43,"PAL") # Pattagansett Lake
text(-72.17, 41.27,"BRI")

dev.off()

