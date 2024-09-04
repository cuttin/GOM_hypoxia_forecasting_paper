if (!require("pacman")) install.packages("pacman"); 
pacman::p_load(magrittr, tidyverse, lubridate, MHadaptive, IDPmisc, 
               reshape2, ggfortify, tictoc, ggthemes, faraway, rstudioapi, gridExtra, ggpubr) 
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
  temp_data <- readRDS(file = paste(path, "ha_new_", yy, ".rds", sep = ""))
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
# Load observations
Out.dat <- read_csv(file = "./data/GoM_cruise_sep18.csv", skip = 1)
Out.dat$Date <- as.Date(Out.dat$Date, format = "%m/%d/%Y")
str(Out.dat)
ha_obs_all <- Out.dat %>% dplyr::select(1, ha_obs = 2, ha_obs_2_5 = LCI, ha_obs_97_5 = UCI,
                          ha_obs_E = 7, ha_sd_E = 9, ha_obs_2_5_E = LCI_1, ha_obs_97_5_E = UCI_1,
                          ha_obs_W = 12, ha_sd_W = 14, ha_obs_2_5_W = LCI_2, ha_obs_97_5_W = UCI_2) 
glimpse(ha_obs_all)
### Load hindcast
hind <- read_csv("./data/hindcast_ha.csv")
glimpse(hind)
pred_dm <- hind %>% 
  dplyr::select(Date,
                ha_pred_e = ha_e_AllUnc_mean,
                ha_pred_w_ba = ha_w_AllUnc_ba_mean) %>% 
  mutate(ha_pred_ba = ha_pred_e + ha_pred_w_ba)
glimpse(pred_dm)
pred_ave <- pred_dm %>%
  mutate(monthday = format(as.Date(Date), "%m-%d")) %>% 
  dplyr::select(-Date) %>% 
  group_by(monthday) %>% 
  summarise_all(mean)
plot(1:122, pred_ave$ha_pred_ba)
pred_dm <- left_join(pred_dm %>% mutate(monthday = format(as.Date(Date), "%m-%d")),
                     pred_ave %>% dplyr::rename(ha_pred_e_ave = ha_pred_e, 
                                                ha_pred_w_ba_ave = ha_pred_w_ba, 
                                                ha_pred_ba_ave = ha_pred_ba),
                     by = 'monthday') %>% 
  dplyr::select(-monthday)
glimpse(pred_dm)
#####################################################################################
## merge all
out_plot <- left_join(outall_m, ha_obs_all, by = "Date") %>% 
  left_join(., pred_dm, by = "Date")
out_plot <- mutate(out_plot, year = lubridate::year(Date))
glimpse(out_plot)
#####################################################################################
### Width of uncertainty, 95% IQR
unc_width <- out_plot %>% 
  mutate(east = `ha_e_AllUnc_97.5%` - `ha_e_AllUnc_2.5%`,
         west = `ha_w_AllUnc_ba_97.5%` - `ha_w_AllUnc_ba_2.5%`) %>% 
  dplyr::select(Date, east, west)
unc_width %>% 
  mutate(ratio = west / east) %>% 
  summarise(mean(ratio))
##### 1 year daily plot
####### By sections, combined
plotit <- out_plot %>% 
  dplyr::select(Date, year, ha_AllUnc_mean = `ha_AllUnc_ba_mean`,
                ha_ParHM_mean = `ha_ParHM_ba_mean`,
                ha_obs = ha_obs,
                ha_obs_2_5 = ha_obs_2_5,
                ha_obs_97_5 = ha_obs_97_5,
                ha_pred = ha_pred_ba, 
                ha_pred_ave = ha_pred_ba_ave,
                ha_AllUnc_2.5 = `ha_AllUnc_ba_2.5%`, ha_AllUnc_97.5 = `ha_AllUnc_ba_97.5%`, 
                ha_ParHMMOD_2.5 = `ha_ParHMMOD_ba_2.5%`, ha_ParHMMOD_97.5 = `ha_ParHMMOD_ba_97.5%`,
                ha_ParHM_2.5 = `ha_ParHM_ba_2.5%`, ha_ParHM_97.5 = `ha_ParHM_ba_97.5%`,
                ha_Par_2.5= `ha_Par_ba_2.5%`, ha_Par_97.5 = `ha_Par_ba_97.5%`) %>% 
  mutate(section = "Combined")
plotit_e <- out_plot %>% 
  dplyr::select(Date, year, ha_AllUnc_mean = `ha_e_AllUnc_mean`,
                ha_ParHM_mean = `ha_e_ParHMUnc_mean`,
                ha_obs = ha_obs_E,
                ha_obs_2_5 = ha_obs_2_5_E,
                ha_obs_97_5 = ha_obs_97_5_E,
                ha_pred = ha_pred_e, 
                ha_pred_ave = ha_pred_e_ave,
                ha_AllUnc_2.5 = `ha_e_AllUnc_2.5%`, ha_AllUnc_97.5 = `ha_e_AllUnc_97.5%`, 
                ha_ParHMMOD_2.5 = `ha_e_ParHMMODUnc_2.5%`, ha_ParHMMOD_97.5 = `ha_e_ParHMMODUnc_97.5%`,
                ha_ParHM_2.5 = `ha_e_ParHMUnc_2.5%`, ha_ParHM_97.5 = `ha_e_ParHMUnc_97.5%`,
                ha_Par_2.5= `ha_e_Par_2.5%`, ha_Par_97.5 = `ha_e_Par_97.5%`) %>% 
  mutate(section = "East")
plotit_w <- out_plot %>% 
  dplyr::select(Date, year, ha_AllUnc_mean = `ha_w_AllUnc_ba_mean`,
                ha_ParHM_mean = `ha_w_ParHMUnc_ba_mean`,
                ha_obs = ha_obs_W,
                ha_obs_2_5 = ha_obs_2_5_W,
                ha_obs_97_5 = ha_obs_97_5_W,
                ha_pred = ha_pred_w_ba,
                ha_pred_ave = ha_pred_w_ba_ave,
                ha_AllUnc_2.5 = `ha_w_AllUnc_ba_2.5%`, ha_AllUnc_97.5 = `ha_w_AllUnc_ba_97.5%`, 
                ha_ParHMMOD_2.5 = `ha_w_ParHMMODUnc_ba_2.5%`, ha_ParHMMOD_97.5 = `ha_w_ParHMMODUnc_ba_97.5%`,
                ha_ParHM_2.5 = `ha_w_ParHMUnc_ba_2.5%`, ha_ParHM_97.5 = `ha_w_ParHMUnc_ba_97.5%`,
                ha_Par_2.5= `ha_w_Par_ba_2.5%`, ha_Par_97.5 = `ha_w_Par_ba_97.5%`) %>% 
  mutate(section = "West")
plotit_all <- bind_rows(plotit, plotit_e) %>% 
  bind_rows(., plotit_w)
##########################################################################################
### Remove negative values
hist(plotit_all %>% dplyr::select(ha_AllUnc_2.5) %>% pull())
## Make negative ones close to zero
plotit_all <- plotit_all %>%
  mutate(ha_AllUnc_2.5 = if_else(ha_AllUnc_2.5 < 0, 0.001, ha_AllUnc_2.5))
##########################################################################################
####################################################### Plot
####
all_unc_ha_plot <- ggplot(plotit_all %>% filter(section == 'Combined', year == 1987))+
  geom_ribbon(aes(x = Date, ymin = ha_AllUnc_2.5, ymax = ha_AllUnc_97.5), fill = "grey30")+
  geom_ribbon(aes(x = Date, ymin = ha_ParHMMOD_2.5, ymax = ha_ParHMMOD_97.5), fill = "grey50")+
  geom_ribbon(aes(x = Date, ymin = ha_ParHM_2.5, ymax = ha_ParHM_97.5), fill = "grey70")+
  geom_ribbon(aes(x = Date, ymin = ha_Par_2.5, ymax = ha_Par_97.5), fill = "grey90")+
  #facet_wrap(~ year, scales = "free_x", ncol = 2)+
  theme_bw()+
  #scale_x_date(expand = c(0, 0), date_labels = "%b",
  #             limits = c(as.Date("2009-06-01"), as.Date("2009-10-01")))+
  scale_y_continuous(limits = c(0,30001), 
                     breaks = c(0, 10000, 20000, 30000), labels = c(0, 10000, 20000, 30000), 
                     expand = c(0, 0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        plot.margin = unit(c(0.5, 0.4, 0, 0.4), "cm"))+
  #geom_line(aes(Date, predicted_w), color = "darkgreen")+
  #geom_line(aes(Date, ha_AllUnc_mean), color = "blue", size = 1.1)+
  #geom_line(aes(Date, ha_pred), color = "#b2abd2", size = 1.1)+
  #geom_line(data = out_plot %>% filter(year == 2013, month == 6), aes(x = Date, y = `ha_AllUnc_mean`), color = "black", size = 1.1)+
  #geom_line(aes(x = Date, y = ha_AllUnc_mean), color = "blue", size = 1.1)+
  #geom_line(aes(x = Date, y = ha_ParHM_mean), color = "red", size = 1.1)+
  geom_point(aes(Date, ha_obs), color = "#e66101", size = 2)+ #
  geom_errorbar(aes(Date, ymin = ha_obs_2_5, ymax = ha_obs_97_5), color = "#e66101", width = 2, size = 0.8)+
  labs(x = "", y = expression(Hypoxic~Area~(km^{2})))
all_unc_ha_plot
########################
########################
plotit_all %>% filter(section != 'Combined') %>%
  group_by(section) %>%
  summarize(max(ha_AllUnc_97.5))
plotit_all$section <- factor(plotit_all$section, 
                                       levels = c("West", "East"), 
                                       ordered = TRUE)
###
##########################################################################################################
#### This one plots two years 1993 and 2009
sel_years <- c(1993, 2009)
for (i in 1:length(sel_years)) {
year_now <- sel_years[i]
ann_text_shelves <- data.frame(Date  = rep(as.Date(paste(year_now, "09", "18", sep = '-')), 2), 
                               yarea =    rep(29000/10000, 2),
                               section = c("West", "East"),
                               label = c(paste0("West, ", year_now), paste0("East, ", year_now)))
ann_text_shelves$section <- factor(ann_text_shelves$section, 
                                   levels = c("West", "East"), 
                                   ordered = TRUE)
plot_one_year_west <- ggplot(plotit_all %>% filter(section == 'West', year == year_now))+
  geom_ribbon(aes(x = Date, ymin = ha_AllUnc_2.5/10000, ymax = ha_AllUnc_97.5/10000), fill = "grey30")+
  geom_ribbon(aes(x = Date, ymin = ha_ParHMMOD_2.5/10000, ymax = ha_ParHMMOD_97.5/10000), fill = "grey50")+
  geom_ribbon(aes(x = Date, ymin = ha_ParHM_2.5/10000, ymax = ha_ParHM_97.5/10000), fill = "grey70")+
  geom_ribbon(aes(x = Date, ymin = ha_Par_2.5/10000, ymax = ha_Par_97.5/10000), fill = "grey90")+
  geom_text(data = ann_text_shelves %>% filter(section == 'West'), 
            aes(x = Date, y = yarea), label = ann_text_shelves$label[ann_text_shelves$section == 'West'], 
            color = "black", size = 3)+
  #facet_wrap(~ section, ncol = 2, scales = "free")+
  theme_bw()+
  #ggtitle(year_now)+
  scale_x_date(expand = c(0.001, 0.001), date_labels = "%b",
               limits = c(as.Date(paste(year_now, "06", "01", sep = '-')), 
                          as.Date(paste(year_now, "10", "01" , sep = "-"))))+
  scale_y_continuous(limits = c(-50/10000,30001/10000), 
                     breaks = c(0, 10000/10000, 20000/10000, 30000/10000), 
                     labels = c(0, 10000/10000, 20000/10000, 30000/10000), 
                     expand = c(0, 0.0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        ## strip 
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(size = 10),
        ## text
        #axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        panel.spacing = unit(1.2, "lines"),
        panel.border = element_rect(colour = "black"),
        #plot.title = element_text(hjust = 0.5, size = 12, margin = margin(0,0,-5,0)),
        #plot.title = element_text(hjust = 0.98, size = 12, vjust = -18), #margin = margin(0,0,0,0)
        plot.margin = unit(c(0.4, 0.3, 0, 0.4), "cm"))+
        #plot.margin = unit(c(-0.5, 0.4, 0, 0.4), "cm"))+
  #geom_line(aes(Date, predicted_w), color = "darkgreen")+
  #geom_line(aes(Date, ha_AllUnc_mean), color = "blue", size = 1.1)+
  geom_line(aes(Date, ha_pred/10000), color = "#fdb863", size = 0.8, linetype = "dashed")+
  geom_line(aes(Date, ha_pred_ave/10000), color = "black", size = 0.8, linetype = "dashed")+
  #geom_line(data = plotit_all %>% 
  #            mutate(month = month(Date)) %>%
  #            filter(section != 'Combined', year == year_now, month == 6), aes(x = Date, y = `ha_AllUnc_mean`), color = "black", size = 1.1, linetype = "dashed")+
  #geom_line(aes(x = Date, y = ha_AllUnc_mean), color = "blue", size = 0.8)+
  #geom_line(aes(x = Date, y = ha_ParHM_mean), color = "red", size = 1.1)+
  geom_point(aes(Date, ha_obs/10000), color = "#e66101", size = 2)+ #
  geom_errorbar(aes(Date, ymin = ha_obs_2_5/10000, ymax = ha_obs_97_5/10000), color = "#e66101", width = 2, size = 0.5)+ #size is width of the vertical bar
  labs(x = "", y = expression(Hypoxic~Area~(10^{4}~km^{2})))
assign(paste('plot_one_year_west', year_now, sep = "_"), plot_one_year_west)
####
#### East
plot_one_year_east <- ggplot(plotit_all %>% filter(section == 'East', year == year_now))+
  geom_ribbon(aes(x = Date, ymin = ha_AllUnc_2.5/10000, ymax = ha_AllUnc_97.5/10000), fill = "grey30")+
  geom_ribbon(aes(x = Date, ymin = ha_ParHMMOD_2.5/10000, ymax = ha_ParHMMOD_97.5/10000), fill = "grey50")+
  geom_ribbon(aes(x = Date, ymin = ha_ParHM_2.5/10000, ymax = ha_ParHM_97.5/10000), fill = "grey70")+
  geom_ribbon(aes(x = Date, ymin = ha_Par_2.5/10000, ymax = ha_Par_97.5/10000), fill = "grey90")+
  geom_text(data = ann_text_shelves %>% filter(section == 'East'), 
            aes(x = Date, y = yarea), label = ann_text_shelves$label[ann_text_shelves$section == 'East'], 
            color = "black", size = 3)+
  #facet_wrap(~ section, ncol = 2, scales = "free")+
  theme_bw()+
  #ggtitle(year_now)+
  scale_x_date(expand = c(0.001, 0.001), date_labels = "%b",
               limits = c(as.Date(paste(year_now, "06", "01", sep = '-')), 
                          as.Date(paste(year_now, "10", "01" , sep = "-"))))+
  scale_y_continuous(limits = c(-50/10000,30001/10000), 
                     breaks = c(0, 10000/10000, 20000/10000, 30000/10000), 
                     labels = c(0, 10000/10000, 20000/10000, 30000/10000), 
                     expand = c(0, 0.0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        ## strip 
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(size = 10),
        ## text
        #axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.ticks.x = element_line(), 
        axis.text.y = element_blank(),
        panel.spacing = unit(1.2, "lines"),
        panel.border = element_rect(colour = "black"),
        #plot.title = element_text(hjust = 0.5, size = 12, margin = margin(0,0,-5,0)),
        #plot.title = element_text(hjust = 0.98, size = 12, vjust = -18), #margin = margin(0,0,0,0)
        plot.margin = unit(c(0.4, 0.4, 0, 0.0), "cm"))+
  #plot.margin = unit(c(-0.5, 0.4, 0, 0.4), "cm"))+
  #geom_line(aes(Date, predicted_w), color = "darkgreen")+
  #geom_line(aes(Date, ha_AllUnc_mean), color = "blue", size = 1.1)+
  geom_line(aes(Date, ha_pred/10000), color = "#fdb863", size = 0.8, linetype = "dashed")+
  geom_line(aes(Date, ha_pred_ave/10000), color = "black", size = 0.8, linetype = "dashed")+
  geom_point(aes(Date, ha_obs/10000), color = "#e66101", size = 2)+ #
  geom_errorbar(aes(Date, ymin = ha_obs_2_5/10000, ymax = ha_obs_97_5/10000), color = "#e66101", width = 2, size = 0.5)+ #size is width of the vertical bar
  labs(x = "", y = "")
assign(paste('plot_one_year_east', year_now, sep = "_"), plot_one_year_east)
}
grid.arrange(plot_one_year_west_1993, plot_one_year_east_1993,
             plot_one_year_west_2009, plot_one_year_east_2009, ncol = 2, widths = c(1.1, 1))

ggsave(filename = paste('./graphics for manuscript/time_series_2_shelfs_combined_1993-2009', '.png', sep = ''),  
       plot = grid.arrange(plot_one_year_west_1993, plot_one_year_east_1993,
                           plot_one_year_west_2009, plot_one_year_east_2009, ncol = 2, widths = c(1.1, 1)), 
       dpi = 500, width = 20,
       height = 16,
       units = c("cm"))
##########################################################################################
##########################################################################################
### Add legend 
legend.data <- data.frame(uncertainty = rep(c("Parameter", "Input data", "Model error", "Adjustment"), 4),
                          x = rep(c(1:4), 4),
                          y = rep(c(1:4), 4)) 
glimpse(legend.data)
legend.data$uncertainty <- factor(legend.data$uncertainty, levels = c("Parameter", "Input data", "Model error", "Adjustment"))
lplot <- ggplot(legend.data)+
  geom_line(aes(x, y, color = uncertainty), size = 5)+
  theme_few()+
  labs(color = "")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.background=element_rect(fill = alpha("white", 0), color = "grey20", size = 0.3),
        legend.key.size = unit(0.3, 'cm'))+
  scale_color_manual(values = c("grey90", "grey70", "grey50", "grey30"))
lplot
legend <- get_legend(lplot)
as_ggplot(legend)
plot_one_year_east_1993_with_legend <- plot_one_year_east_1993+
  annotation_custom(legend, xmin = as.Date(paste(1993, "08", "15", sep = '-')), ymin = 1.2)
plot_one_year_east_1993_with_legend 
grid.arrange(plot_one_year_west_1993, plot_one_year_east_1993_with_legend,
             plot_one_year_west_2009, plot_one_year_east_2009, ncol = 2, widths = c(1.1, 1))
##########################################################################################
##########################################################################################
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## This plot with the legend is the one in the manuscipt
ggsave(filename = paste('./graphics for manuscript/time_series_2_shelfs_combined_1993-2009_with_legend', '.png', sep = ''),  
       plot = grid.arrange(plot_one_year_west_1993, plot_one_year_east_1993_with_legend,
                           plot_one_year_west_2009, plot_one_year_east_2009, ncol = 2, widths = c(1.1, 1)), 
       dpi = 500, width = 20,
       height = 16,
       units = c("cm"))
##########################################################################################
##########################################################################################
##########################################################################################
#### Plot all years by shelf 
years_of_calibration <- 1985:2016
plotlist_west <- list()
plotlist_east <- list()
for (i in seq_along(years_of_calibration)) {
  
  yearsoi <- years_of_calibration[i]
  ann_text_shelves <- data.frame(Date  = rep(as.Date(paste(yearsoi, "09", "18", sep = '-')), 2), 
                                 yarea =    rep(29000/10000, 2),
                                 section = c("West", "East"),
                                 label = c(paste0("West, ", yearsoi), paste0("East, ", yearsoi)))
  ann_text_shelves$section <- factor(ann_text_shelves$section, 
                                     levels = c("West", "East"), 
                                     ordered = TRUE)
  plot_one_year_west <- ggplot(plotit_all %>% filter(section == 'West', year == yearsoi))+
    geom_ribbon(aes(x = Date, ymin = ha_AllUnc_2.5/10000, ymax = ha_AllUnc_97.5/10000), fill = "grey30")+
    geom_ribbon(aes(x = Date, ymin = ha_ParHMMOD_2.5/10000, ymax = ha_ParHMMOD_97.5/10000), fill = "grey50")+
    geom_ribbon(aes(x = Date, ymin = ha_ParHM_2.5/10000, ymax = ha_ParHM_97.5/10000), fill = "grey70")+
    geom_ribbon(aes(x = Date, ymin = ha_Par_2.5/10000, ymax = ha_Par_97.5/10000), fill = "grey90")+
    geom_text(data = ann_text_shelves %>% filter(section == 'West'), 
              aes(x = Date, y = yarea), label = ann_text_shelves$label[ann_text_shelves$section == 'West'], 
              color = "black", size = 3)+
    #facet_wrap(~ section, ncol = 2, scales = "free")+
    theme_bw()+
    #ggtitle(yearsoi)+
    scale_x_date(expand = c(0.001, 0.001), date_labels = "%b",
                 limits = c(as.Date(paste(yearsoi, "06", "01", sep = '-')), 
                            as.Date(paste(yearsoi, "10", "01" , sep = "-"))))+
    scale_y_continuous(limits = c(-50/10000,30001/10000), 
                       breaks = c(0, 10000/10000, 20000/10000, 30000/10000), 
                       labels = c(0, 10000/10000, 20000/10000, 30000/10000), 
                       expand = c(0, 0.0))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          ## strip 
          strip.text.x = element_blank(),
          strip.background = element_blank(),
          #strip.text = element_text(size = 10),
          ## text
          #axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          panel.spacing = unit(1.2, "lines"),
          panel.border = element_rect(colour = "black"),
          #plot.title = element_text(hjust = 0.5, size = 12, margin = margin(0,0,-5,0)),
          #plot.title = element_text(hjust = 0.98, size = 12, vjust = -18), #margin = margin(0,0,0,0)
          plot.margin = unit(c(0.4, 0.3, 0, 0.4), "cm"))+
    geom_line(aes(Date, ha_pred/10000), color = "#fdb863", size = 0.8, linetype = "dashed")+
    geom_line(aes(Date, ha_pred_ave/10000), color = "black", size = 0.8, linetype = "dashed")+
    geom_point(aes(Date, ha_obs/10000), color = "#e66101", size = 2)+ #
    geom_errorbar(aes(Date, ymin = ha_obs_2_5/10000, ymax = ha_obs_97_5/10000), color = "#e66101", width = 2, size = 0.5)+ #size is width of the vertical bar
    labs(x = "", y = expression(Hypoxic~Area~(10^{4}~km^{2})))
  plot_one_year_west
  ####
  #### East
  plot_one_year_east <- ggplot(plotit_all %>% filter(section == 'East', year == yearsoi))+
    geom_ribbon(aes(x = Date, ymin = ha_AllUnc_2.5/10000, ymax = ha_AllUnc_97.5/10000), fill = "grey30")+
    geom_ribbon(aes(x = Date, ymin = ha_ParHMMOD_2.5/10000, ymax = ha_ParHMMOD_97.5/10000), fill = "grey50")+
    geom_ribbon(aes(x = Date, ymin = ha_ParHM_2.5/10000, ymax = ha_ParHM_97.5/10000), fill = "grey70")+
    geom_ribbon(aes(x = Date, ymin = ha_Par_2.5/10000, ymax = ha_Par_97.5/10000), fill = "grey90")+
    geom_text(data = ann_text_shelves %>% filter(section == 'East'), 
              aes(x = Date, y = yarea), label = ann_text_shelves$label[ann_text_shelves$section == 'East'], 
              color = "black", size = 3)+
    #facet_wrap(~ section, ncol = 2, scales = "free")+
    theme_bw()+
    #ggtitle(yearsoi)+
    scale_x_date(expand = c(0.001, 0.001), date_labels = "%b",
                 limits = c(as.Date(paste(yearsoi, "06", "01", sep = '-')), 
                            as.Date(paste(yearsoi, "10", "01" , sep = "-"))))+
    scale_y_continuous(limits = c(-50/10000,30001/10000), 
                       breaks = c(0, 10000/10000, 20000/10000, 30000/10000), 
                       labels = c(0, 10000/10000, 20000/10000, 30000/10000), 
                       expand = c(0, 0.0))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          ## strip 
          strip.text.x = element_blank(),
          strip.background = element_blank(),
          #strip.text = element_text(size = 10),
          ## text
          #axis.title = element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.ticks.x = element_line(), 
          axis.text.y = element_blank(),
          panel.spacing = unit(1.2, "lines"),
          panel.border = element_rect(colour = "black"),
          #plot.title = element_text(hjust = 0.5, size = 12, margin = margin(0,0,-5,0)),
          #plot.title = element_text(hjust = 0.98, size = 12, vjust = -18), #margin = margin(0,0,0,0)
          plot.margin = unit(c(0.4, 0.4, 0, 0.0), "cm"))+
    geom_line(aes(Date, ha_pred/10000), color = "#fdb863", size = 0.8, linetype = "dashed")+
    geom_line(aes(Date, ha_pred_ave/10000), color = "black", size = 0.8, linetype = "dashed")+
    geom_point(aes(Date, ha_obs/10000), color = "#e66101", size = 2)+ #
    geom_errorbar(aes(Date, ymin = ha_obs_2_5/10000, ymax = ha_obs_97_5/10000), color = "#e66101", width = 2, size = 0.5)+ #size is width of the vertical bar
    labs(x = "", y = "")
  plot_one_year_east
  
  plotlist_west[[i]] <- plot_one_year_west
  plotlist_east[[i]] <- plot_one_year_east
}


date_now <- Sys.Date()
pdf(paste("./graphics for manuscript/HA_by_shelf_", date_now, ".pdf", sep = ""), 
    onefile = TRUE, width = 8, height = 5)
for (i in 1:length(plotlist_west)) {
  grid.arrange(plotlist_west[[i]], plotlist_east[[i]], ncol = 2, widths = c(1.1, 1))
}
dev.off()

## Save png
date_now <- Sys.Date()
iter <- seq(1, 31, by = 2)
iter <- seq(3, 32, by = 3)
iter <- c(1, iter)
for (i in iter) {
  if (i < 2) {
  ggsave(filename = paste('./graphics for manuscript/HA_shelf_', date_now, "_", i, ".png", sep = ""),  
         
         plot = grid.arrange(plotlist_west[[i]], plotlist_east[[i]], 
                             plotlist_west[[i+1]], plotlist_east[[i+1]], ncol = 2, widths = c(1.1, 1)), 
         dpi = 500, width = 20,
         height = 15,
         units = c("cm")) } else {
  
           ggsave(filename = paste('./graphics for manuscript/HA_shelf_', date_now, "_", i, ".png", sep = ""),  
                  
                  plot = grid.arrange(plotlist_west[[i]], plotlist_east[[i]], 
                                      plotlist_west[[i+1]], plotlist_east[[i+1]], 
                                      plotlist_west[[i+2]], plotlist_east[[i+2]], ncol = 2, widths = c(1.1, 1)), 
                  dpi = 500, width = 20,
                  height = 22.5,
                  units = c("cm"))         
                    
         }
}
dev.off()
#########################################################################################################################################
#########################################################################################################################################
### Combined area
years_of_calibration <- 1985:2016
plotlist <- list()
every4 <- seq(1, 32, by = 4)

for (i in seq_along(every4)) {
  ii = every4[i]
  yearsoi <- years_of_calibration[ii:(ii+3)]
  ann_text_years <- data.frame(Date  = c(as.Date(paste(yearsoi[1], "09", "23", sep = '-')),
                                         as.Date(paste(yearsoi[2], "09", "23", sep = '-')),
                                         as.Date(paste(yearsoi[3], "09", "23", sep = '-')),
                                         as.Date(paste(yearsoi[4], "09", "23", sep = '-'))), 
                              yarea = 39000 / 10000,
                               year =  yearsoi,
                              label = yearsoi)
  
  p <- ggplot(plotit %>% filter(year %in% yearsoi))+
    #ggtitle(yearsoi)+
    geom_ribbon(aes(x = Date, ymin = ha_AllUnc_2.5/10000, ymax = ha_AllUnc_97.5/10000), fill = "grey30")+
    geom_ribbon(aes(x = Date, ymin = ha_ParHMMOD_2.5/10000, ymax = ha_ParHMMOD_97.5/10000), fill = "grey50")+
    geom_ribbon(aes(x = Date, ymin = ha_ParHM_2.5/10000, ymax = ha_ParHM_97.5/10000), fill = "grey70")+
    geom_ribbon(aes(x = Date, ymin = ha_Par_2.5/10000, ymax = ha_Par_97.5/10000), fill = "grey90")+
    facet_wrap(~ year, scales = "free", ncol = 2)+
    theme_bw()+
    geom_text(data = ann_text_years, 
              aes(x = Date, y = yarea), label = ann_text_years$label, 
              color = "black", size = 4)+
    #scale_x_date(expand = c(0.001, 0.001), date_labels = "%b",
    #             limits = c(as.Date(paste(yearsoi, "06", "01", sep = '-')), 
    #                        as.Date(paste(yearsoi, "10", "01" , sep = "-"))))+
    scale_y_continuous(limits = c(-50/10000, 40000/10000), 
                       breaks = c(0, 10000/10000, 20000/10000, 30000/10000), 
                       labels = c(0, 10000/10000, 20000/10000, 30000/10000), 
                       expand = c(0, 0.0))+
    theme_bw()+
    scale_x_date(expand = c(0, 0), date_labels = "%b")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(), # removes strip text
          panel.spacing = unit(1.2, "lines"),
          panel.border = element_rect(colour = "black"),
          plot.margin = unit(c(0.5, 0.4, 0.3, 0.4), "cm"))+
    # theme(panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       strip.background = element_blank(),
    #       strip.text = element_text(size = 10),
    #       axis.title = element_text(size = 10),
    #       axis.text = element_text(size = 8),
    #       panel.spacing = unit(1.2, "lines"),
    #       panel.border = element_rect(colour = "black"),
    #       plot.title = element_text(size = 12),
    #       #plot.title = element_text(hjust = 0.5, size = 12, margin = margin(0,0,-5,0)),
    #       #plot.title = element_text(hjust = 0.98, size = 12, vjust = -18), #margin = margin(0,0,0,0)
    #       plot.margin = unit(c(-0.5, 0.4, 0, 0.4), "cm"))+
    #geom_line(aes(Date, predicted_w), color = "darkgreen")+
    #geom_line(aes(Date, ha_AllUnc_mean), color = "blue", size = 1.1)+
    geom_line(aes(Date, ha_pred/10000), color = "#fdb863", size = 0.8, linetype = "dashed")+
    geom_line(aes(Date, ha_pred_ave/10000), color = "black", size = 0.8, linetype = "dashed")+
    #geom_line(data = plotit %>% 
    #            mutate(month = month(Date)) %>%
    #            filter(section != 'Combined', year == year_now, month == 6), aes(x = Date, y = `ha_AllUnc_mean`), color = "black", size = 1.1, linetype = "dashed")+
    #geom_line(aes(x = Date, y = ha_AllUnc_mean), color = "blue", size = 0.8)+
    #geom_line(aes(x = Date, y = ha_ParHM_mean), color = "red", size = 1.1)+
    geom_point(aes(Date, ha_obs/10000), color = "#e66101", size = 2)+ #
    geom_errorbar(aes(Date, ymin = ha_obs_2_5/10000, ymax = ha_obs_97_5/10000), color = "#e66101", width = 2, size = 0.5)+ #size is width of the vertical bar
    labs(x = "", y = expression(Hypoxic~Area~(10^{4}~km^{2})))
  p
  plotlist[[i]] <- p
}



date_now <- Sys.Date()
pdf(paste("./graphics for manuscript/HA_total_", date_now, ".pdf", sep = ""), 
    onefile = TRUE, width = 8, height = 9)
for (i in 1:length(plotlist)) {
  print(plotlist[[i]])
}
dev.off()

# save png 
date_now <- Sys.Date()

for (i in 1:length(plotlist)) {
  ggsave(filename = paste('./graphics for manuscript/HA_total_', date_now, "_", i, ".png", sep = ""),  
         plot = plotlist[[i]], 
         dpi = 500, width = 20,
         height = 20,
         units = c("cm"))
  #print(plotlist[[i]])
}
dev.off()
####
####
#### Observations outside of predictive interval
####
out_in <- data.frame(Date = plotit$Date)
out_in$outside <- if_else(plotit$ha_obs > plotit$ha_AllUnc_97.5 | plotit$ha_obs < plotit$ha_AllUnc_2.5, 1, 0)
length(na.omit(out_in$outside))
sum(out_in$outside, na.rm = T)
# percent outside the predictive interval
sum(out_in$outside, na.rm = T) / length(na.omit(out_in$outside)) * 100
na.omit(out_in$Date[out_in$outside == 1])
####
####
#### Uncertainty bands difference
glimpse(plotit)
unc <- plotit %>% 
  mutate(All = ha_AllUnc_97.5 - ha_AllUnc_2.5,
         ParHM = ha_ParHM_97.5 - ha_ParHM_2.5,
         Par = ha_Par_97.5 - ha_Par_2.5,
         HM = ParHM - Par,
         ErrPar = All - HM,
         Err = All - ParHM, # predictive error and due to  transformation
         month = month(Date)) %>%
  dplyr::select(Date, month, year, All, ParHM, ErrPar, Par, HM, Err)
unc$monthday <- format(as.Date(unc$Date), "%m-%d")
unc$time <- if_else(unc$month == 6 | unc$month == 7, "early", "late")
summary(unc$time)
which(unc$time == "early")
ggplot(unc, aes(monthday, HM, color = month))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(unc, aes(monthday, (HM)/Par, color = month))+
  geom_boxplot()
ggplot(unc, aes(monthday, (HM)/Err, color = month))+
  geom_boxplot()

# Riverine-meteo and both parameter-error
unc %>% 
  mutate(parerr_hm_diff = HM / ErrPar) %>%
  summarise(mean(parerr_hm_diff))
unc %>% 
  mutate(parerr_hm_diff = HM / ErrPar) %>%
  group_by(time) %>%
  summarise(mean(parerr_hm_diff))
# Riverine-meteo and parameter
unc %>% 
  mutate(par_hm_diff = HM / Par) %>%
  summarise(mean(par_hm_diff))
unc %>% 
  mutate(par_hm_diff = HM / Par) %>%
  group_by(time) %>%
  summarise(mean(par_hm_diff))

# Riverine-meteo and error
unc %>% 
  mutate(err_hm_diff = HM / Err) %>%
  summarise(mean(err_hm_diff))
unc %>% 
  mutate(err_hm_diff = HM / Err) %>%
  group_by(time) %>%
  summarise(mean(err_hm_diff))

# Parameter and error uncertainty
unc %>% 
  mutate(par_err_diff = Err / Par) %>%
  summarise(mean(par_err_diff))
unc %>% 
  mutate(par_err_diff = Err / Par) %>%
  group_by(time) %>%
  summarise(mean(par_err_diff))


unc %>% 
  mutate(par_hm_diff = HM / Par) %>%
  group_by(time) %>%
  summarise(median(par_hm_diff))
unc %>% 
  mutate(err_hm_diff = (HM - Err) / Err) %>%
  group_by(month) %>%
  summarise(mean(err_hm_diff))

unc %>% 
  mutate(par_hm_diff = (HM - Par) / Par) %>%
  group_by(month) %>%
  summarise(mean(par_hm_diff))
unc %>% 
  mutate(err_hm_diff = (HM - Err) / Err) %>%
  group_by(month) %>%
  summarise(mean(err_hm_diff))
