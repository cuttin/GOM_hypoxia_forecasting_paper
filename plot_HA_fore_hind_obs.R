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
years_of_interest <- c(1985:2016)
for (i in seq_along(years_of_interest)) {
  yy <- years_of_interest[i]
  temp_data <- readRDS(file = paste(path, "ha_", yy, ".rds", sep = ""))
  ha_data <- bind_rows(ha_data, temp_data)
}
glimpse(ha_data)
## Extract from this one hm, parameter uncertainty
ha_data_old <- data.frame()
years_of_interest <- c(1985:2016)
for (i in seq_along(years_of_interest)) {
  yy <- years_of_interest[i]
  temp_data <- readRDS(file = paste(path, "ha_old_", yy, ".rds", sep = ""))
  ha_data_old <- bind_rows(ha_data_old, temp_data)
}
glimpse(ha_data_old)
####
outall_m <- bind_cols(ha_data_old %>% 
                        dplyr::select(-contains("AllUnc")),
                      ha_data %>% 
                        dplyr::select(contains("AllUnc")))
dplyr::glimpse(outall_m)
#####################################################################################
## Get hindcast values
hind <- read_csv("./data/hindcast_ha.csv")
glimpse(hind)
pred_dm <- hind %>% 
  dplyr::select(Date,
                ha_pred_e = ha_e_AllUnc_mean,
                ha_pred_w_ba = ha_w_AllUnc_ba_mean) %>% 
  mutate(ha_pred_ba = ha_pred_e + ha_pred_w_ba)
glimpse(pred_dm)
#####################################################################################
## Get observed values
Out.dat <- read_csv(file = "./data/GoM_cruise_sep18.csv", skip = 1)
Out.dat$Date <- as.Date(Out.dat$Date, format = "%m/%d/%Y")
str(Out.dat)
ha_obs_all <- Out.dat %>% dplyr::select(1, ha_obs = 2, ha_obs_2_5 = LCI, ha_obs_97_5 = UCI,
                                        ha_obs_E = 7, ha_sd_E = 9, ha_obs_2_5_E = LCI_1, ha_obs_97_5_E = UCI_1,
                                        ha_obs_W = 12, ha_sd_W = 14, ha_obs_2_5_W = LCI_2, ha_obs_97_5_W = UCI_2) 
glimpse(ha_obs_all)
#####################################################################################
## Combined forec vs hind and obs
out_plot_fho <- left_join(outall_m %>% dplyr::select(Date, Forecasted = ha_AllUnc_ba_mean),
                          pred_dm  %>% dplyr::select(Date, Hindcasted = ha_pred_ba), by = 'Date') %>% 
  left_join(., ha_obs_all %>% dplyr::select(Date, Observed = ha_obs), by = 'Date') %>% 
  mutate(month = month(Date))
out_plot_fho_m <- melt(out_plot_fho, id.vars = c("Date", "month", "Forecasted"))
labels_fho <- c(expression(Hindcasted~(10^{4}~km^{2})),
                expression(Observed~(10^{4}~km^{2})))
glimpse(out_plot_fho_m)
out_plot_fho_m$variable <- factor(out_plot_fho_m$variable, 
                                  levels = c("Hindcasted", "Observed"), 
                                  ordered = TRUE, labels = labels_fho)
hindcasted <- as.character(unique(out_plot_fho_m$variable)[1])
observed <- as.character(unique(out_plot_fho_m$variable)[2])
ann_text_r2_fho <- data.frame(value = c(5000, 5000),
                              Forecasted = c(22000, 22000),
                              variable = c(hindcasted, observed),
                              label = c(paste("italic(R) ^ 2 == ", 
                                              round(E(out_plot_fho_m$Forecasted[out_plot_fho_m$variable == hindcasted], 
                                                      out_plot_fho_m$value[out_plot_fho_m$variable == hindcasted]), 2), 
                                              sep = ""),
                                        paste("italic(R) ^ 2 == ", 
                                              round(E(out_plot_fho_m$Forecasted[out_plot_fho_m$variable == observed], 
                                                      out_plot_fho_m$value[out_plot_fho_m$variable == observed]), 2), 
                                              sep = "")))
ggplot(out_plot_fho_m , aes(y = value, x = Forecasted))+
  geom_point(shape = 21)+
  facet_wrap( ~ variable , strip.position = "left", labeller = label_parsed)+
  geom_abline(slope = 1)+
  theme_few()+
  geom_text(data = ann_text_r2_fho, label = ann_text_r2_fho$label, parse = T, color = "black")+
  labs(x = expression(Forecasted~(km^{2})), y = "")+
  theme(strip.placement = "outside",
        strip.text = element_text(size = 10), # text size inside strip (y-axis title)
        axis.title.x = element_text(size = 10), # text size inside strip (x-axis title) 
        axis.text = element_text(size = 8)) # text size axis values
###############################
out_plot_fho_month_m <- out_plot_fho_m
out_plot_fho_month_m$monthf <- month.abb[out_plot_fho_month_m$month]
out_plot_fho_month_m$monthf <- parse_factor((out_plot_fho_month_m$monthf), levels = month.abb)
(hfha6 <- format(round(E(out_plot_fho_month_m$Forecasted[out_plot_fho_month_m$variable == hindcasted & out_plot_fho_month_m$month == 6], 
                         out_plot_fho_month_m$value[out_plot_fho_month_m$variable == hindcasted & out_plot_fho_month_m$month == 6]), 2), nsmall = 2))
(hfha7 <- format(round(E(out_plot_fho_month_m$Forecasted[out_plot_fho_month_m$variable == hindcasted & out_plot_fho_month_m$month == 7], 
                         out_plot_fho_month_m$value[out_plot_fho_month_m$variable == hindcasted & out_plot_fho_month_m$month == 7]), 2), nsmall = 2))
(hfha8 <- format(round(E(out_plot_fho_month_m$Forecasted[out_plot_fho_month_m$variable == hindcasted & out_plot_fho_month_m$month == 8], 
                         out_plot_fho_month_m$value[out_plot_fho_month_m$variable == hindcasted & out_plot_fho_month_m$month == 8]), 2), nsmall = 2))
(hfha9 <- format(round(E(out_plot_fho_month_m$Forecasted[out_plot_fho_month_m$variable == hindcasted & out_plot_fho_month_m$month == 9], 
                         out_plot_fho_month_m$value[out_plot_fho_month_m$variable == hindcasted & out_plot_fho_month_m$month == 9]), 2), nsmall = 2))  
(ofha6 <- format(round(E(out_plot_fho_month_m$Forecasted[out_plot_fho_month_m$variable == observed & out_plot_fho_month_m$month == 6], 
                         out_plot_fho_month_m$value[out_plot_fho_month_m$variable == observed & out_plot_fho_month_m$month == 6]), 2), nsmall = 2))
(ofha7 <- format(round(E(out_plot_fho_month_m$Forecasted[out_plot_fho_month_m$variable == observed & out_plot_fho_month_m$month == 7], 
                         out_plot_fho_month_m$value[out_plot_fho_month_m$variable == observed & out_plot_fho_month_m$month == 7]), 2), nsmall = 2))
(ofha8 <- format(round(E(out_plot_fho_month_m$Forecasted[out_plot_fho_month_m$variable == observed & out_plot_fho_month_m$month == 8], 
                         out_plot_fho_month_m$value[out_plot_fho_month_m$variable == observed & out_plot_fho_month_m$month == 8]), 2), nsmall = 2))
(ofha9 <- format(round(E(out_plot_fho_month_m$Forecasted[out_plot_fho_month_m$variable == observed & out_plot_fho_month_m$month == 9], 
                         out_plot_fho_month_m$value[out_plot_fho_month_m$variable == observed & out_plot_fho_month_m$month == 9]), 2), nsmall = 2))
## Combined forecasted versus hindcast and observed, by month
ann_text_r2_fho_month <- data.frame(value = rep(2000, 4),
                                    Forecasted = rep(25000, 4),
                                    variable = c(rep(hindcasted, 4),
                                                 rep(observed, 4)),
                                    monthf = rep(unique(out_plot_fho_month_m$monthf), 2),
                                    # label = c(paste("italic(R) ^ 2 == ", hfha6, sep = ""),
                                    #           paste("italic(R) ^ 2 == ", hfha7, sep = ""),
                                    #           paste("italic(R) ^ 2 == ", hfha8, sep = ""),
                                    #           paste("italic(R) ^ 2 == ", hfha9, sep = ""),
                                    #           
                                    #           paste("italic(R) ^ 2 == ", ofha6, sep = ""),
                                    #           paste("italic(R) ^ 2 == ", ofha7, sep = ""),
                                    #           paste("italic(R) ^ 2 == ", ofha8, sep = ""),
                                    #           paste("italic(R) ^ 2 == ", ofha9, sep = "")))
                                    label = c(hfha6, hfha7, hfha8, hfha9, ofha6, ofha7, ofha8, ofha9))

ann_text_r2_fho_month$monthf <- month.abb[ann_text_r2_fho_month$monthf]
ann_text_r2_fho_month$monthf <- parse_factor(ann_text_r2_fho_month$monthf, levels = month.abb)
###
monthly_ha_fho <- ggplot(out_plot_fho_month_m, aes(y = value/10000, x = Forecasted/10000))+ #, color = monthf
  geom_point(shape = 21)+
  facet_rep_grid(variable ~ monthf, labeller = label_parsed, switch = "y", repeat.tick.labels = F)+
  geom_abline(slope = 1)+
  theme_few()+
  geom_text(data = ann_text_r2_fho_month[ann_text_r2_fho_month$variable == hindcasted & ann_text_r2_fho_month$monthf == "Jun",], 
            label = paste0("italic(R) ^ 2 == ", deparse(ann_text_r2_fho_month$label[1])), parse = T, color = "black", size = 3)+
  geom_text(data = ann_text_r2_fho_month %>% filter(variable != hindcasted | monthf != "Jun"), 
            label = paste0("italic(R) ^ 2 == ", ann_text_r2_fho_month$label[2:8]), parse = T, color = "black", size = 3)+
  #geom_text(data = ann_text_r2_fho_month, label = ann_text_r2_fho_month$label, parse = T, color = "black", size = 3)+
  labs(x = expression(Forecasted~HA~(10^{4}~km^{2})), y = "", color = "")+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 31000/10000), 
                     #breaks = c(0, 5000, 15000, 25000))+
                     breaks = c(0, 10000/10000, 20000/10000, 30000/10000))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 31000/10000), 
                     #breaks = c(0, 5000, 15000, 25000))+
                     breaks = c(0, 10000/10000, 20000/10000, 30000/10000))+
  #scale_color_brewer(palette = "Spectral")+
  theme(legend.position = "none",
    strip.placement = "outside",
        strip.text = element_text(size = 10), # text size inside strip (y-axis title)
        axis.title.x = element_text(size = 10), # text size inside strip (x-axis title) 
        axis.text = element_text(size = 8)) # text size axis values
monthly_ha_fho
time_now <- Sys.Date()
png(filename = paste("./graphics for manuscript/hof_monthly_HA_ba_", time_now, ".png", sep = ""), 
    units = "cm", 
    width = 20, 
    height = 12, 
    #pointsize=12, 
    res = 500)
print(monthly_ha_fho)
dev.off()
###############################################################################################
# Just obs, hind vs forecasted r-squared
E(out_plot_fho$Forecasted, out_plot_fho$Hindcasted)
E(out_plot_fho$Forecasted, out_plot_fho$Observed)
###############################################################################################
## Combined forecasted versus hindcast and observed, by shelf
out_plot_fhoE <- left_join(outall_m %>% dplyr::select(Date, Forecasted = ha_e_AllUnc_mean), 
                           pred_dm %>% dplyr::select(Date, Hindcasted = ha_pred_e), by = "Date") %>% 
              left_join(., ha_obs_all %>% dplyr::select(Date, Observed = ha_obs_E), by = "Date") %>% 
  mutate(shelf = "East")
out_plot_fhoW <- left_join(outall_m %>% dplyr::select(Date, Forecasted = ha_w_AllUnc_ba_mean), 
                           pred_dm %>% dplyr::select(Date, Hindcasted = ha_pred_w_ba), by = "Date") %>% 
              left_join(., ha_obs_all %>% dplyr::select(Date, Observed = ha_obs_W), by = "Date") %>% 
  mutate(shelf = "West")
out_plot_fho_shelves <- bind_rows(out_plot_fhoE, out_plot_fhoW)
#######
#######
glimpse(out_plot_fho_shelves)
#######
#######
out_plot_fho_shelves_m <- melt(out_plot_fho_shelves, id.vars = c("Date", "shelf", "Forecasted"))

out_plot_fho_shelves_m$variable <- factor(out_plot_fho_shelves_m$variable, 
                                          levels = c("Hindcasted", "Observed"), 
                                          ordered = TRUE, labels = labels_fho)
out_plot_fho_shelves_m$shelf <- factor(out_plot_fho_shelves_m$shelf, 
                                       levels = c("West", "East"), 
                                       ordered = TRUE)
ann_text_r2_fho_shelves <- data.frame(value = rep(2000, 4),
                                      Forecasted = rep(20000, 4),
                                      variable = c(hindcasted,
                                                   hindcasted,
                                                   observed,
                                                   observed),
                                      shelf = c("West", "East", "West", "East"),
                                      label = c(paste("italic(R) ^ 2 == ", 
                                                      round(E(out_plot_fhoW$Forecasted, out_plot_fhoW$Hindcasted), 2), sep = ""),
                                                paste("italic(R) ^ 2 == ", 
                                                      round(E(out_plot_fhoE$Forecasted, out_plot_fhoE$Hindcasted), 2), sep = ""),
                                                paste("italic(R) ^ 2 == ", 
                                                      round(E(out_plot_fhoW$Forecasted, out_plot_fhoW$Observed), 2), sep = ""),
                                                paste("italic(R) ^ 2 == ", 
                                                      round(E(out_plot_fhoE$Forecasted, out_plot_fhoE$Observed), 2), sep = "")))
ann_text_r2_fho_shelves$shelf <- factor(ann_text_r2_fho_shelves$shelf, 
                                        levels = c("West", "East"), 
                                        ordered = TRUE)
out_plot_fho_shelves_m <- mutate(out_plot_fho_shelves_m, month = month(Date))
out_plot_fho_shelves_m$monthf <- month.abb[out_plot_fho_shelves_m$month]
out_plot_fho_shelves_m$monthf <- parse_factor((out_plot_fho_shelves_m$monthf), levels = month.abb)
hfo_ha <- ggplot(out_plot_fho_shelves_m, aes(y = value/10000, x = Forecasted/10000, color = monthf))+ #, color = monthf
  geom_point(shape = 21)+
  facet_rep_grid(variable ~ shelf, switch = "y", labeller = label_parsed, repeat.tick.labels = F)+
  geom_abline(slope = 1)+
  theme_few()+
  geom_text(data = ann_text_r2_fho_shelves, label = ann_text_r2_fho_shelves$label, 
            parse = T, color = "black", size = 3)+
  #scale_color_manual(values = Palette_months())+
  scale_color_brewer(palette = "Spectral")+
  labs(x = expression(Forecasted~(10^{4}~km^{2})), y = "", color = "")+
  scale_x_continuous(expand = c(0.05, 0.05), limits = c(0, 25000/10000), breaks = c(0, 10000/10000, 20000/10000))+
  scale_y_continuous(expand = c(0.05, 0.05), limits = c(0, 25000/10000), breaks = c(0, 10000/10000, 20000/10000))+
  theme(legend.position = "top",
        legend.text = element_text(size = 8),
        strip.placement = "outside",
        strip.text = element_text(size = 10), # text size inside strip (y-axis title)
        axis.title.x = element_text(size = 10), # text size inside strip (x-axis title) 
        axis.text = element_text(size = 8), # text size axis values
        legend.margin = margin(0,0,0,0),
        panel.spacing = unit(1, "lines"),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"), #top,right,bot, left
        #axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.9),
        legend.box.margin = margin(-7,15,-5,-7))
hfo_ha
time_now <- Sys.Date()
png(filename = paste("./graphics for manuscript/hind_fore_two_shelves_", time_now, ".png", sep = ""),
    units = "cm", 
    width = 14, 
    height = 12,
    res = 500)
print(hfo_ha)
dev.off()
#########################################################