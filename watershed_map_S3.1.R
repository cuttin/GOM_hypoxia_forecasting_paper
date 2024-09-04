if (!require("pacman")) install.packages("pacman"); 
pacman::p_load(adaptMCMC, TTR, truncnorm, tidyverse, lubridate, suncalc, data.table, deSolve,
               zoo, faraway, gridExtra, reshape2, DataCombine, leaps, gdata, ggthemes, 
               cowplot, scales, GGally, odin, gtable, rstudioapi, matrixStats, hydroGOF, EGRET,
               rnoaa, sp, raster, gstat, animation, ggmap, tmap, sf, ggspatial, rgdal, ggsn, png,
               leaflet, prism, devtools, maps, mapproj, usmap)
select <- dplyr::select
basepath <- "Q:\\My Drive\\research\\transfered from Katin\\GOM model\\Data"
setwd(basepath)
getwd()
#####################################################################################
#####################################################################################
# from here https://www.sciencebase.gov/catalog/item/55de04d5e4b0518e354dfcf8
# Read watershed coordinates 
#dem <- raster::raster("Q:\\My Drive\\research\\transfered from Katin\\GOM model\\Data\\miswatershed\\c0_PolygonToRaster1.tif")
dem <- raster::raster("Q:\\My Drive\\research\\transfered from Katin\\GOM model\\Data\\misswatershedfromsciencebase\\Miss_RiverBasin\\Miss_RiverBasin_PolygonToRas1.tif")
plot(dem)
sr <- "+proj=longlat"
demp <- projectRaster(dem, crs = sr)
demdf <- as(demp, 'SpatialGridDataFrame')
df <- as.data.frame(demdf)
dff <- df %>% select(longitude = s1, latitude = s2)
str(dff)
############# River
mars <- st_read("Q:\\My Drive\\research\\transfered from Katin\\GOM model\\Data\\MajorRiver_ERF1\\Major_River\\Major_River.shp")
mars <- st_crs(mars, crs = sr)
marsraster <- rasterize(mars, demp)
#############
### Map of states
us_states <- map_data("state")
us_states$states <- str_to_title(us_states$region)
us_states$abbreviation <- state.abb[match(us_states$states, state.name)]
head(us_states$abbreviation)
centroid <- aggregate(cbind(long,lat) ~ abbreviation, data=us_states, FUN=mean)
centroid.xy <- us_states %>% 
  group_by(abbreviation) %>%
  summarise(long = mean(range(long)), 
            lat = mean(range(lat)),
            label = unique(abbreviation))



plot_map_1 <- ggplot()+ 
  geom_polygon(data = us_states,
               mapping = aes(x = long, y = lat,
                             group = group), color = "grey20", size = 0.1, fill = "white") +
  geom_point(data = dff, aes(x = longitude, y = latitude, group = 1), color = 'lightblue', alpha = 0.05)+
  #geom_text(data = centroid, mapping = aes(x=long, y=lat, label=abbreviation, group = 1), size = 2)+
  geom_text(data = centroid.xy, mapping = aes(x=long, y=lat, label=label, group = 1), size = 2)+
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  guides(fill = FALSE)+
  theme_map()

plot_map_1

png("plot_map_1.png", res = 500, width = 10, height = 10, units = "cm")
print(plot_map_1)
dev.off()

plot_usmap(labels = TRUE)
