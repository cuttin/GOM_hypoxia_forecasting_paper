if (!require("pacman")) install.packages("pacman"); 
pacman::p_load(adaptMCMC, TTR, truncnorm, tidyverse, lubridate, suncalc, data.table, deSolve,
               zoo, faraway, gridExtra, reshape2, DataCombine, leaps, gdata, ggthemes, 
               cowplot, scales, GGally, odin, gtable, rstudioapi, matrixStats, hydroGOF, EGRET,
               rnoaa, sp, raster, gstat, animation, ggmap, tmap, sf, ggspatial, rgdal, ggsn, png,
               leaflet, prism, devtools, maps, mapproj, usmap, akima, metR, ggpattern)
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
centroid.xy <- us_states %>% 
  group_by(abbreviation) %>%
  summarise(long = mean(range(long)), 
            lat = mean(range(lat)),
            label = unique(abbreviation))
centroid.xy$long[centroid.xy$label == "FL"] <- -81.2
centroid.xy$long[centroid.xy$label == "LA"] <- -92.4
centroid.xy$long[centroid.xy$label == "MI"] <- -85.1
centroid.xy$long[centroid.xy$label == "ID"] <- -114.5

# ### This one is used to print the names of the states
# centroid.xy <- bind_rows(
#   bind_cols(long = -86.7, lat = 32.6, label = "AL"),
#  # bind_cols(long = -81.2, lat = 28.1, label = "FL"),
#   bind_cols(long = -92.4, lat = 31.0,  label = "LA"),
#   bind_cols(long = -89.8, lat = 32.6, label = "MS"),
#   bind_cols(long = -99.0, lat = 31.2, label = "TX"))
### Read seamap data to get the polygon
### Get hypoxia
jun_matli <- readr::read_csv('Q:\\My Drive\\research\\transfered from Katin\\!GoM shrimps\\DO_2020M\\JuneDO.csv')
### Extract one day data to create polygon
glimpse(jun_matli)
grid <- jun_matli %>% 
  filter(date == unique(jun_matli$date)[100]) %>%
  select(east, north) 
gridplot <- grid
coordinates(grid) <- ~ east + north
CRSstring <- paste0("+proj=utm +zone=", "15N ", "+units=km")
utmcoor <- sp::SpatialPoints(grid, proj4string = sp::CRS(CRSstring))
utmcoordf <- data.frame(utmcoor)
longlatcoor <- sp::spTransform(utmcoor, sp::CRS("+init=epsg:4326"))
lonlatdf <- bind_cols(lon = data.frame(longlatcoor)$east, 
                      lat = data.frame(longlatcoor)$north)
ggplot(gridplot, aes(east, north))+geom_point()
#################################################################################################
(latext <- c(min(lonlatdf$lat) + 0.1, 31)) #max(lonlatdf$lat) + 1.5))
(lonext <- c(max(lonlatdf$lon) + 0.2, min(lonlatdf$lon) - 0.2))
baseUrldark<- "https://api.mapbox.com/styles/v1/vmatli/cjwuu3bco80mn1cqe9m6lb4vf/tiles/256/{z}/{x}/{y}?access_token=pk.eyJ1Ijoidm1hdGxpIiwiYSI6ImNqczR5ZmMxeTBhNngzeXJ3bGZ0czJuZ3QifQ.VSDH4msUykRloXReWVu_-A"
##  mapbox://styles/vmatli/cjwjmem0x00aw1cnq6in8v2m5  mapbox://styles/vmatli/cjwuu3bco80mn1cqe9m6lb4vf
##"https://api.mapbox.com/styles/v1/mapbox/dark-v10/tiles/256/{z}/{x}/{y}?access_token=pk.eyJ1Ijoidm1hdGxpIiwiYSI6ImNqczR5ZmMxeTBhNngzeXJ3bGZ0czJuZ3QifQ.VSDH4msUykRloXReWVu_-A"
myMaposmd <- OpenStreetMap::openmap(upperLeft = c(latext[2], lonext[2]),
                                    lowerRight = c(latext[1], lonext[1]),
                                    type = baseUrldark,zoom = 9.1)
sr <- "+proj=longlat"
myMaposmd <- OpenStreetMap::openproj(myMaposmd, projection = sr)
#################################################################################################
# Read Mis atch rivers
misatch <- sf::st_read('Q:\\My Drive\\research\\transfered from Katin\\GOM model\\Clean\\data\\study site map\\MissAtch.shp')
ggplot() + 
  geom_sf(data = misatch, size = 1, color = "black") + 
  ggtitle("Mis & Atch") + 
  coord_sf()
str(misatch)
(misatch_df <- as.data.frame(st_coordinates(st_transform(misatch, crs = sr))))
misatch_df <- bind_cols(as.data.frame((st_transform(misatch, crs = sr))),
                        as.data.frame(st_coordinates(st_transform(misatch, crs = sr))))
hist(misatch_df$L1)
glimpse(misatch_df)
Atch <- data.frame(x = -91.7, y = 30.6, label = "Atchafalaya River")
Miss <- data.frame(x = -90.1, y = 30.6, label = "Mississippi River")
#################################################################################################
# Read depth contours
contours <- sf::st_read('Q:\\My Drive\\research\\transfered from Katin\\GOM model\\Clean\\data\\study site map\\Contours10SelectSmooth2.shp')
ggplot() + 
  geom_sf(data = contours, size = 1, color = "black") + 
  ggtitle("Mis & Atch") + 
  coord_sf()
str(contours)
contours_df <- as.data.frame(st_coordinates(st_transform((contours), crs = sr))) 
glimpse(contours_df)
hist(contours_df$L1)
## Depth
depth <- sf::st_read('Q:\\My Drive\\research\\transfered from Katin\\GOM model\\Clean\\data\\study site map\\EstPointDepth2.shp')
ggplot() + 
  geom_sf(data = depth, size = 1, color = "black") + 
  ggtitle("Mis & Atch") + 
  coord_sf()
str(depth)
unique(depth$ShelfSect)
depth_df <- as.data.frame(st_coordinates(st_transform(depth, crs = sr))) 
depth_df <- bind_cols(as.data.frame(depth), 
                      as.data.frame(st_coordinates(st_transform(depth, crs = sr))) %>%
                        rename(lon = X, lat = Y)) %>% 
  drop_na()
glimpse(depth_df) 
pnts_west <- depth_df %>%
  filter(ShelfSect == 1) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
polygon_west <- concaveman::concaveman(pnts_west, concavity = 10) #, concavity = 4
plot(polygon_west, reset = FALSE)
poly_west <- data.frame(st_coordinates(polygon_west[,1])) %>% dplyr::select(lon = X, lat = Y) %>% 
  mutate(section = "West")
pnts_east <- depth_df %>%
  filter(ShelfSect == 2) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
polygon_east <- concaveman::concaveman(pnts_east, concavity = 10)
poly_east <- data.frame(st_coordinates(polygon_east[,1])) %>% dplyr::select(lon = X, lat = Y) %>% 
  mutate(section = "East")
poly <- bind_rows(poly_west, poly_east)
glimpse(poly)
#################################################################################################
# Grid from depth
grid <- with(depth_df, interp(x = lon, y = lat, z = RASTERVALU, linear = TRUE, extrap = TRUE,
                              xo = seq(min(lon), max(lon), length = 100), 
                              yo = seq(min(lat), max(lat), length = 100))) 
ggplot(depth_df, aes(lon, lat, color = RASTERVALU))+geom_point(size = 4)
glimpse(depth_df)
dfll <- as.data.frame(interp2xyz(grid))
ggplot(dfll, aes(x=x, y=y))+ 
  geom_contour(aes(z = z), breaks = c(-20, -30, -40, -50, -60, -70), linetype = "dashed", color = "grey50")+
  metR::geom_text_contour(aes(z = z), rotate = FALSE, breaks = c(-20, -30, -40, -50, -60, -70), label.placer = label_placer_flattest())
#################################################################################################
## Read arrow
arrow <- readPNG("Q:\\My Drive\\research\\transfered from Katin\\Data\\north-arrow.png")
g_arrow <- rasterGrob(arrow, interpolate = TRUE)
#################################################################################################
########################################################################################
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Finally plot
gommap <- OpenStreetMap::autoplot.OpenStreetMap(myMaposmd, expand=FALSE)+
  #geom_point(data = lonlatdf, aes(lon, lat))+
  geom_point(data = misatch_df %>% filter(Y < 31, Y > 29.5), aes(X, Y, group = L1), color = cbbPalette[c(6)], size = 0.7)+
  geom_polygon(data = poly, aes(lon, lat, fill = section, color = section), size = 2.5)+ # 
  theme_bw()+
  labs(fill = '', color = '', x = '', y = '')+ #x = 'lon', y = "lat"
  scale_y_continuous(expand = c(0.00, 0.00), limits = c(min(lonlatdf$lat) - 0.05, 30.7))+
  scale_x_continuous(expand = c(0.03, 0.03))+
  geom_text(data = Miss, aes(x = x, y = y, label = label), size = 4, color = "black", #cbbPalette[c(6)], 
            hjust = 1.2, nudge_x = 0.003)+
  geom_text(data = Atch, aes(x = x, y = y, label = label), size = 4, color = "black", #cbbPalette[c(6)],
            hjust = 1.2, nudge_x = 0.003)+
  #geom_line(data = contours_df, aes(X,Y, group = L1))+
  geom_contour(data = dfll, aes(x = x, y = y, z = z), 
               breaks = c(-20, -30, -40, -50, -60, -70), linetype = "dashed", color = "grey50")+
  metR::geom_text_contour(data = dfll, aes(x = x, y = y, z = z), 
                          rotate = FALSE, breaks = c(-20, -30, -40, -50, -60, -70), 
                          label.placer = label_placer_flattest(),
                          size = 3)+
  scale_fill_manual(values = cbbPalette[c(2, 4)])+
  scale_color_manual(values = cbbPalette[c(2, 4)])+
  theme(legend.position = "top", #c(0, 1)
        legend.justification = c(0, 1),
        strip.background = element_blank(),
        #legend.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(1,1,1,1),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 12),
        #text = element_text(size = 8),
        #legend.title = element_text(size = 12),
        plot.margin = unit(c(-0.1,0.1,0.1,0.1), "cm"))+
  annotation_custom(g_arrow, xmin = (lonext[1]-0.3), xmax = (lonext[1]), 
                    ymin = (latext[2]-0.2), ymax = (latext[2]-0.05))+
  scalebar(x.min = (lonext[2]), x.max = (lonext[1]-0.18),
           y.min = (latext[1]-0.03), y.max = (latext[2]),
           #height = .2, 
           st.size = 2,
           transform = T,
           dist = 25, dist_unit = "km", model = "WGS84") 
  
gommap
########################################################################################
########################################################################################
########################################################################################
us <- ggplot()+ 
  geom_polygon(data = us_states,
               mapping = aes(x = long, y = lat,
                             group = group), color = "grey20", size = 0.1, fill = "white") +
  geom_point(data = dff, aes(x = longitude, y = latitude, group = 1), color = 'lightblue', alpha = 0.05)+
  #geom_text(data = centroid.xy %>% filter(!(label %in% c('NE', 'MA', 'RI', 'CT', 'NJ', 'DE', 'MD'))), mapping = aes(x=long, y=lat, label=label, group = 1), size = 2)+
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  theme_minimal()+
  geom_rect(mapping = aes(xmin = lonext[1]+1, xmax = lonext[2], 
                          ymin = latext[1]-0.3, ymax = latext[2]-0.4), 
            color = cbbPalette[c(7)], fill = NA, size = 1)+
  #geom_segment(aes(x = (lonext[1]+1+lonext[2])/2, y = latext[1]-0.3, 
  #                 xend = (lonext[1]+1+lonext[2])/2, yend = latext[1]-4.5),
  #             lineend = "butt", linejoin = "round", color = cbbPalette[c(7)],
  #             size = 1.2, arrow = arrow(length = unit(0.1, "inches")))+
  theme(plot.margin = unit(c(-0.1,-0.5,-0.1,-0.5), "cm"), #t = 0, r = 0, b = 0, l = 0
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(colour = "black", fill="white", size=0.5))
us
########################################################################################
########################################################################################
########################################################################################
finplot <- gommap + annotation_custom(ggplotGrob(us), xmin = lonext[2]+0.05, xmax = lonext[2]+1.5, 
                           ymin = latext[2]-1.4, ymax = latext[2]-0.1)
finplot
#ggpubr::ggarrange(us, gommap, nrow = 2, ncol = 1)

(date_now <- Sys.Date())
ggsave(filename = paste0('Q:\\My Drive\\research\\transfered from Katin\\GOM model\\Clean\\data\\study site map\\fig_1_v2_', date_now, '.pdf'),  
       plot = finplot,   
       dpi = 500, width = 20, height = 10, units = c("cm"))
#######################################################################################
