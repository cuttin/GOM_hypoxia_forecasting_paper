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
#############
### Map of states
us_states <- map_data("state")
us_states$states <- str_to_title(us_states$region)
us_states$abbreviation <- state.abb[match(us_states$states, state.name)]
glimpse(us_states)
southern_us_states <- filter(us_states, abbreviation %in% c("TX", "LA", "MS", "AL", "FL"))
southern_us_states <- filter(us_states, abbreviation %in% c("TX", "LA", "MS", "AL"))
head(us_states$abbreviation)
centroid <- aggregate(cbind(long,lat) ~ abbreviation, data=southern_us_states, FUN=mean)
centroid.xy <- southern_us_states %>% 
  group_by(abbreviation) %>%
  summarise(long = mean(range(long)), 
            lat = mean(range(lat)),
            label = unique(abbreviation))
centroid.xy

plot_usmap(include = .south_region, labels = TRUE)
plot_usmap(include = c("TX", "LA", "MS", "AL", "FL"), labels = TRUE)

### This one is used to print the names of the states
centroid.xy <- bind_rows(
  bind_cols(long = -86.7, lat = 32.6, label = "AL"),
 # bind_cols(long = -81.2, lat = 28.1, label = "FL"),
  bind_cols(long = -92.4, lat = 31.0,  label = "LA"),
  bind_cols(long = -89.8, lat = 32.6, label = "MS"),
  bind_cols(long = -99.0, lat = 31.2, label = "TX"))
### Read seamap data to get the polygon
df_all <- read.table('Q:\\My Drive\\research\\transfered from Katin\\!GoM shrimps\\SEAMAP\\Shrimp codes\\brown shrimp\\bshrimp binomial/sum.lu.brown.csv', sep=',', header =TRUE)
df_all$Date <- as.Date(paste(df_all$YR, df_all$month, df_all$day), "%Y %m %d")
df_all <- df_all %>% arrange(Date)
glimpse(df_all)
#################################################################################################
(latext <- c(min(dat$lat) + 0.1, max(dat$lat) + 0.8))
(lonext <- c(max(dat$lon) + 0.1, min(dat$lon) - 0.2))
baseUrldark<- "https://api.mapbox.com/styles/v1/vmatli/cjwuu3bco80mn1cqe9m6lb4vf/tiles/256/{z}/{x}/{y}?access_token=pk.eyJ1Ijoidm1hdGxpIiwiYSI6ImNqczR5ZmMxeTBhNngzeXJ3bGZ0czJuZ3QifQ.VSDH4msUykRloXReWVu_-A"
##  mapbox://styles/vmatli/cjwjmem0x00aw1cnq6in8v2m5  mapbox://styles/vmatli/cjwuu3bco80mn1cqe9m6lb4vf
##"https://api.mapbox.com/styles/v1/mapbox/dark-v10/tiles/256/{z}/{x}/{y}?access_token=pk.eyJ1Ijoidm1hdGxpIiwiYSI6ImNqczR5ZmMxeTBhNngzeXJ3bGZ0czJuZ3QifQ.VSDH4msUykRloXReWVu_-A"
myMaposmd <- OpenStreetMap::openmap(upperLeft = c(latext[2], lonext[2]),
                                    lowerRight = c(latext[1], lonext[1]),
                                    type = baseUrldark,zoom = 9.1)
sr <- "+proj=longlat"
myMaposmd <- OpenStreetMap::openproj(myMaposmd, projection = sr)
#################################################################################################
## Read arrow
arrow <- readPNG("Q:\\My Drive\\research\\transfered from Katin\\Data\\north-arrow.png")
g_arrow <- rasterGrob(arrow, interpolate = TRUE)
#################################################################################################
########################################################################################
##Create a polygon file from data to limit NASA
pnts_west <- df_all %>%
  filter(lon < -91.2) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
polygon_west <- concaveman::concaveman(pnts_west, concavity = 7)
plot(polygon_west, reset = FALSE)
plot(pnts, add = TRUE)
class(polygon_west)
polygon_west
poly_west <- data.frame(st_coordinates(polygon_west[,1])) %>% dplyr::select(lon = X, lat = Y) %>% 
  mutate(section = "West")
##
pnts_east <- df_all %>%
  filter(lon >= -91.2, lon <= -90) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
polygon_east <- concaveman::concaveman(pnts_east, concavity = 7)
plot(polygon_east, reset = FALSE)
plot(pnts, add = TRUE)
class(polygon_east)
polygon_east
poly_east <- data.frame(st_coordinates(polygon_east[,1])) %>% dplyr::select(lon = X, lat = Y) %>% 
  mutate(section = "East")
poly <- bind_rows(poly_west, poly_east)
##

###
plot_map_1 <- ggplot()+ 
  geom_polygon(data = southern_us_states,
               mapping = aes(x = long, y = lat,
                             group = group), color = "grey20", size = 0.1, fill = "grey80") +
  #geom_point(data = dff, aes(x = longitude, y = latitude, group = 1), color = 'lightblue', alpha = 0.05)+
  #geom_text(data = centroid, mapping = aes(x=long, y=lat, label=abbreviation, group = 1), size = 2)+
  geom_text(data = centroid.xy, mapping = aes(x=long, y=lat, label=label, group = 1), size = 2)+
  geom_polygon(data = poly, aes(lon, lat, fill = section))+
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)
  
  #guides(fill = FALSE)+
  #theme_map()

plot_map_1
#######################################################################################



#######################################################################################
png("plot_map_1.png", res = 500, width = 10, height = 10, units = "cm")
print(plot_map_1)
dev.off()

plot_usmap(labels = TRUE)
