if (!require("pacman")) install.packages("pacman"); 
pacman::p_load(magrittr, tidyverse, lubridate, MHadaptive, IDPmisc, 
               reshape2, ggfortify, tictoc, ggthemes, faraway, rstudioapi) 
select <- dplyr::select
#### Get directory - location of the current R file and set it as basepath
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))
getwd()
#####################################################################################
# All functions here
source("all_functions.R")
#####################################################################################
## Load forecast
## Extract from this one only bias adjustment and transformation uncertainty
path <- paste0(getwd(), "/data/all_sim_09-09_paral/")
ha_data <- data.frame()
ha_data_old <- data.frame()
ha_hm <- data.frame()
years_of_interest <- c(1985:2016)
for (i in seq_along(years_of_interest)) {
  yy <- years_of_interest[i]
  temp_data <- readRDS(file = paste(path, "ha_new_", yy, ".rds", sep = ""))
  ha_data <- bind_rows(ha_data, temp_data)
  temp_data_old <- readRDS(file = paste(path, "ha_old_", yy, ".rds", sep = ""))
  ha_data_old <- bind_rows(ha_data_old, temp_data_old)
  temp_data_hm <- readRDS(file = paste(path, "hm_old_", yy, ".rds", sep = ""))
  ha_hm <- bind_rows(ha_hm, temp_data_hm)
}
glimpse(ha_data)
glimpse(ha_data_old)
glimpse(ha_hm)
###################################################################################################
ggplot(ha_hm %>% mutate(year = year(Date)))+
  geom_ribbon(aes(x = Date, ymin = `ha_hm_ba_2.5%`, ymax = `ha_hm_ba_97.5%`))+
  facet_wrap(~ year, scales = "free_x")
ggplot(ha_hm %>% mutate(monthday = format(as.Date(Date), "%m-%d")))+
  geom_boxplot(aes(x = monthday, y = `ha_hm_ba_97.5%`-`ha_hm_ba_2.5%`))
###########################################################################################################
## Compare with hindcast uncertainty
hind_unc <- read_csv("./data/hindcast_var.csv")
glimpse(hind_unc)
cc <- bind_cols(ha_data %>% dplyr::select(Date, Forecast = ha_AllUnc_ba_var),
                hind_unc %>% dplyr::select(Hindcast = `Parameter, Residual error, Transformation`))
mean(cc$Forecast) / (mean(ha_hm$ha_hm_ba_var) + mean(cc$Hindcast))
cc %>% mutate(month = month(Date)) %>%
ggplot(., aes(Forecast/10e6, Hindcast/10e6, color = factor(month)))+
  geom_point()+
  #facet_wrap(~month)+
  #scale_x_continuous(labels = c(0.1, 0.2, 0.3), breaks = c(0.1, 0.2, 0.3))+
  #scale_y_continuous(labels = c(0.1, 0.2, 0.3), breaks = c(0.1, 0.2, 0.3))+
  coord_fixed(ratio = 1, xlim = c(0, 3.5), ylim = c(0, 3.5), expand = TRUE)+
  ggtitle('Variance comparison')+
  scale_color_brewer(palette = "Spectral")+
  labs(y = expression(Hindcasted~(km^{4})~10^{6}), x = expression(Forecasted~(km^{4})~10^{6}), color = "Month")+
  theme(legend.position=c(.08,.85))+
  geom_abline(slope = 1)
###########################################################################################################
####
ha <- left_join(ha_data_old %>%
                  dplyr::select(Date, Parameter = ha_Par_ba_var, ha_ParHM_ba_var), 
                ha_data %>%
                  dplyr::select(Date, `All uncertainty` = ha_AllUnc_ba_var)) %>%
  dplyr::mutate(`Residual error` = `All uncertainty` - ha_ParHM_ba_var) %>%
  dplyr::select(Date, Parameter, `Residual error`, `All uncertainty`) %>%
  left_join(., ha_hm   %>% dplyr::select(Date, `Data inputs` = ha_hm_ba_var),
            by = 'Date')


unc <- ha %>% mutate(month = month(Date))
unc$time <- if_else(unc$month == 6 | unc$month == 7, "early", "late")
# Riverine-meteo and parameter
unc %>% 
  mutate(err_hm_diff = `Data inputs` / Parameter) %>%
  summarise(mean(err_hm_diff))
unc %>% 
  mutate(err_hm_diff = `Data inputs` / Parameter) %>%
  group_by(time) %>%
  summarise(mean(err_hm_diff))
# Riverine-meteo and residual-transformation
unc %>% 
  mutate(err_hm_diff = `Data inputs` / `Residual error`) %>%
  summarise(mean(err_hm_diff))
unc %>% 
  mutate(err_hm_diff = `Residual error` / `Data inputs`) %>%
  summarise(median(err_hm_diff))
unc %>% 
  summarise(mean(`Residual error`) / mean(`Data inputs`))
unc %>% 
  mutate(err_hm_diff = `Data inputs` / `Residual error`) %>%
  group_by(time) %>%
  summarise(mean(err_hm_diff))

unc %>% 
  mutate(err_hm_diff = `Residual error` /`Data inputs` ) %>%
  group_by(time) %>%
  summarise(mean(err_hm_diff))


# Riverine-meteo and residual-transformation+parameter
unc %>% 
  mutate(err_hm_diff = `Data inputs` / (`Residual error`+Parameter)) %>%
  group_by(time) %>%
  summarise(mean(err_hm_diff))

# Riverine-meteo and all uncertainty
unc %>% 
  mutate(err_hm_diff = `Data inputs` / (`All uncertainty`)) %>%
  group_by(time) %>%
  summarise(mean(err_hm_diff))

##
ham <- ha %>% dplyr::select(-`All uncertainty`) %>%
  melt(., id.vars = "Date")
unique(ham$variable)
glimpse(ham)
ham$monthday <- format(as.Date(ham$Date), "%m-%d")
ggplot(ham %>% filter(variable != "residual error"), aes(monthday, value, color = variable))+
  geom_boxplot(outlier.shape = NA)+
  theme_few()+
  theme(legend.position = "top")+
  #scale_color_brewer(palette = 'Spectral')+
  #scale_fill_brewer(palette = 1)+
  labs(x = "", y = "Variance", color = "")+
  scale_x_discrete(breaks = c("06-01", "06-15", 
                              "07-01", "07-15", 
                              "08-01", "08-15", 
                              "09-01", "09-15"))
## Relative size of variances
## Mean variance
ham$variable <- factor(ham$variable, levels = c("Parameter", 
                                                "Data inputs",
                                                "Residual error"))
var_plot <- ham %>% 
  dplyr::select(-Date) %>%
  group_by(monthday, variable) %>%
  summarise(value = (mean(value))) %>%
  ggplot(., aes(monthday, value/10e6, color = variable, group = variable))+
  geom_line(size = 1.1)+
  theme_few()+
  scale_color_manual(values = Palette_variances())+
  theme(legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(0, 0, 0, 0),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8))+
        scale_y_continuous(expand = c(0.001, 0.001), breaks = c(0, 0.5, 1, 1.5), 
                           labels = c(0, 0.5, 1, 1.5), limits = c(0, 1.5))+
  labs(x = "", y = expression(Averaged~daily~variance~(km^{4}~10^{6})), color = "")+
  scale_x_discrete(breaks = c("06-15", 
                              "07-01", "07-15", 
                              "08-01", "08-15", 
                              "09-01", "09-15"), 
                   labels = c("Jun 15", 
                              "Jul 1", "Jul 15", 
                              "Aug 1", "Aug 15", 
                              "Sep 1", "Sep 15"))
var_plot
date_now <- Sys.Date()
png(filename = paste0("./graphics for manuscript/variance_", date_now, ".png"), 
    units="cm", 
    width = 12, #14
    height = 6, #8.5
    #pointsize=12, 
    res = 300)
print(var_plot)
dev.off()
########################################################################################################