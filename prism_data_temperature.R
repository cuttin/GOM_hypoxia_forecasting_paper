if (!require("pacman")) install.packages("pacman"); 
pacman::p_load(adaptMCMC, TTR, truncnorm, tidyverse, lubridate, suncalc, data.table, deSolve,
               zoo, faraway, gridExtra, reshape2, DataCombine, leaps, gdata, ggthemes, rstan, 
               cowplot, scales, GGally, odin, gtable, rstudioapi, matrixStats, hydroGOF, EGRET,
               rnoaa, sp, raster, gstat, animation, ggmap, tmap, sf, ggspatial, rgdal, ggsn, png,
               leaflet, prism, devtools)
#### Get directory - location of the current R file and set it as basepath
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))
getwd()
#####################################################################################
#####################################################################################
# from here https://www.sciencebase.gov/catalog/item/55de04d5e4b0518e354dfcf8
# Read watershed coordinates 
dem <- raster::raster("./data/MAR watershed raster/Miss_RiverBasin_PolygonToRas1.tif")
plot(dem)
sr <- "+proj=longlat"
demp <- projectRaster(dem, crs = sr)
demdf <- as(demp, 'SpatialGridDataFrame')
df <- as.data.frame(demdf)
dff <- df %>% select(longitude = s1, latitude = s2)
str(dff)
#####################################################################################
#prism::prism_set_dl_dir(/prismppt")
#options(prism.path = "C:\\Users\\akatin\\prismppt")
#options(prism.path = "/prismppt")
options(prism.path = "C:/Temp/prismtemp/")
get_prism_monthlys(type = "tmean", years = c(1980:2016), mon = c(1:12), keepZip = T) # download data 
ls_prism_data(name = TRUE) # see the list of stats
temp <- ls_prism_data(name = TRUE)
temp$product_name
# create a dataframe with year-month from 1985 to 2016
monthly_data <- data.frame(Date = seq.Date(from = as.Date("1980-01-01"), 
                                           to = as.Date("2016-12-31"), by = "month")) %>%
  mutate(year = year(Date), month = month(Date)) %>% 
  select(- Date)
monthly_data$t <- NA
## Get mean monthly temperature
for (i in seq_along(temp$files)) {
  RS <- prism_stack(ls_prism_data()[i, 1]) # save 1 year one month
  proj4string(RS) <- proj4string(demdf) # similar projection as watershed
  RSS <- as(RS, 'SpatialPointsDataFrame') # convert to Spatial Points
  inout <- over(RSS, demdf, returnList = F) %>% pull(Miss_RiverBasin_PolygonToRas1) # Extract the ones inside watershed
  df <- data.frame(rasterToPoints(RS)) ##creates a dataframe of points
  colnames(df) <- c("x", "y", "t")
  df$inout <- inout
  df$inout <- if_else(is.na(df$inout) == TRUE, 0, df$inout)
  d <- df %>% filter(inout == 1) # leave only points within watershed
  monthly_data$t[i] <- mean(d$t) # write to dataframe
  print(i)
}

glimpse(monthly_data)
write_csv(monthly_data, "./data/prism_temp_monthly.csv") ## all data

summary(monthly_data$t)
win <- monthly_data %>% 
  filter(month %in% c(1:3)) %>%
  group_by(year) %>%
  summarise(jfm_T = mean(t, na.rm = T))  
spr <- monthly_data %>% 
  filter(month %in% c(4:5)) %>%
  group_by(year) %>%
  summarise(am_T = mean(t, na.rm = T))
may <- monthly_data %>% 
  filter(month %in% c(5)) %>%
  group_by(year) %>%
  summarise(m_T = mean(t, na.rm = T))
temperature <- left_join(win, spr, by = "year") %>% 
  left_join(., may, by = "year")
summary(temperature)
write_csv(temperature, "./data/prism_temp_monthly_watershed.csv")
###################
###################
###################
###################
################### Steps of process (don't run)
###################
###################
new_file <- c(1) ##change to corresponding file numbers
RS <- prism_stack(ls_prism_data()[new_file, 1])
proj4string(RS) <- proj4string(demdf)
RSS <- as(RS, 'SpatialPointsDataFrame') # convert to Spatial Points
inout <- over(RSS, demdf, returnList = F) %>% pull(c0_PolygonToRaster1) # Extract the ones inside watershed
df <- data.frame(rasterToPoints(RS)) ##creates a dataframe of points
str(df)
colnames(df) <- c("x", "y", "ppt")
df$inout <- inout
df$inout <- if_else(is.na(df$inout) == TRUE, 0, df$inout)
dff <- df %>% filter(inout == 1)
str(dff)
mean(dff$ppt)
#
plot1 <- ggplot(data = dff, aes(x = longitude, y = latitude))+
  geom_point(color = 'grey80')+
  geom_point(data = df, aes(x = x, y = y, color = factor(inout)), alpha = 0.3)+
  theme(legend.position = "bottom")+
  labs(color = "Station within watershed: 0-out, 1-in", x = "", y = "")
png("checkitt.png", height = 8, width = 10, units = "cm")
plot1
dev.off()
