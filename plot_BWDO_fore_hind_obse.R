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
path <- paste0(getwd(), "/data/all_sim_09-09_paral/")
fore_bwdo <- data.frame()
years_of_interest <- c(1985:2016)
for (i in seq_along(years_of_interest)) {
  yy <- years_of_interest[i]
  temp_data <- readRDS(file = paste(path, "bwdo_", yy, ".rds", sep = ""))
  fore_bwdo <- bind_rows(fore_bwdo, temp_data)
}
dplyr::glimpse(fore_bwdo)
#####################################################################################
## Load hindcast 
hind <- read_csv("./data/hindcast_bwdo.csv")
glimpse(hind)
#####################################################################################
## Get observed values
Out.dat <- read_csv(file = "./data/GoM_cruise_sep18.csv", skip = 1)
Out.dat$Date <- as.Date(Out.dat$Date, format = "%m/%d/%Y")
glimpse(Out.dat)
Out.dat.all <- Out.dat %>% 
  dplyr::select(Date, obsOw = Mean_4,   obsOe = Mean_3, # define observations 
                errOw = STDEV_4,  errOe = STDEV_3) # define sdv geostat errors 
glimpse(Out.dat.all)
#####################################################################################
## Combine
ploto <- bind_rows(left_join(fore_bwdo %>% dplyr::select(Date, Forecasted = bwdoE_e_mean),
                                  hind %>% dplyr::select(Date, Hindcasted = bwdoE_hind)  , by = 'Date') %>%
              left_join(., Out.dat.all %>% dplyr::select(Date, Observed   = obsOe)       , by = 'Date') %>% 
              mutate(shelf = "East"),
              ## West shelf
              left_join(fore_bwdo %>% dplyr::select(Date, Forecasted = bwdoW_e_ba_mean),
                        hind %>% dplyr::select(Date, Hindcasted = bwdoW_hind)  , by = 'Date') %>%
                left_join(., Out.dat.all %>% dplyr::select(Date, Observed   = obsOw)       , by = 'Date') %>%
                     mutate(shelf = "West")) %>% 
  mutate(month = month(Date), year = year(Date))
ploto$shelf <- factor(ploto$shelf, levels = c("West", "East"), labels = c("West", "East"))
summary(ploto$month)
ploto$monthf <- month.abb[ploto$month]
ploto$monthf <- parse_factor((ploto$monthf), levels = month.abb)
glimpse(ploto)
#####################################################################################
## Obs vs pred
(r2e_op <- round(E(ploto$Hindcasted[ploto$shelf == "East"], ploto$Observed[ploto$shelf == "East"]), 2))
(r2w_op <- round(E(ploto$Hindcasted[ploto$shelf == "West"], ploto$Observed[ploto$shelf == "West"]), 2))
(r2_op <- round(E(ploto$Hindcasted, ploto$Observed), 2))
## Pred forec
(r2e_pf <- round(E(ploto$Forecasted[ploto$shelf == "East"], ploto$Hindcasted[ploto$shelf == "East"]), 2))
(r2w_pf <- round(E(ploto$Forecasted[ploto$shelf == "West"], ploto$Hindcasted[ploto$shelf == "West"]), 2))
(r2_pf <- round(E(ploto$Forecasted, ploto$Hindcasted), 2))
## Obs forec
(r2e_of <- round(E(ploto$Forecasted[ploto$shelf == "East"], ploto$Observed[ploto$shelf == "East"]), 2))
(r2w_of <- round(E(ploto$Forecasted[ploto$shelf == "West"], ploto$Observed[ploto$shelf == "West"]), 2))
(r2_of <- round(E(ploto$Forecasted, ploto$Observed), 2))
#####################################################################################
ann_text_pf_ew <- data.frame(Forecasted = c(4.5), Hindcasted = c(2.5),
                             shelf = c('East', 'West'),
                             label = c(paste("italic(R) ^ 2 == ", r2e_pf, sep = ""),
                                       paste("italic(R) ^ 2 == ", r2w_pf, sep = "")))
## Hindcasted - Forec
ann_text_pf <- data.frame(Forecasted = c(4.5), Hindcasted = c(2.5),
                          label = c(paste("italic(R) ^ 2 == ", r2_pf, sep = "")))
ggplot(ploto, aes(y = Hindcasted, Forecasted, color = shelf))+
  geom_point(shape = 21, alpha = 0.7)+
  #facet_wrap(~ shelf, ncol = 2)+
  geom_abline(slope = 1)+
  theme_few()+
  theme(legend.position = "top")+
  labs(color = "")+
  geom_text(data = ann_text_pf, label = ann_text_pf$label, parse = T, color = "black")
## Hind - forec by shelf
ann_text_pfs <- data.frame(Forecasted = c(5.2, 5.2), Hindcasted = c(1.5, 1.5),
                           shelf = c("West", "East"),
                           label = c(paste("italic(R) ^ 2 == ", r2w_pf, sep = ""),
                                     paste("italic(R) ^ 2 == ", r2e_pf, sep = "")))
ann_text_pfs$shelf <- factor(ann_text_pfs$shelf, levels = c("West", "East"))
hind_fore_plot_two_shelves <- ggplot(ploto, aes(y = Hindcasted, Forecasted, color = monthf))+
  geom_point(shape = 21, alpha = 1)+
  facet_grid(.~ shelf)+
  geom_abline(slope = 1, color = "grey50", size = 0.3)+
  theme_few()+
  theme(legend.position = "top")+
  scale_color_manual(values = Palette_months())+
  coord_fixed(ratio = 1, xlim = c(1, 6), ylim = c(1, 6), expand = TRUE)+
  labs(color = "", x = expression(Forecasted~BWDO~(mg/L)), y = expression(Hindcasted~BWDO~(mg/L)))+
  geom_text(data = ann_text_pfs, label = ann_text_pfs$label, parse = T, color = "black", size = 3)+
  theme(strip.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,-10,-10))
hind_fore_plot_two_shelves
####################################################################################
#### Monthly performance degradation
ann_text_pfsm <- data.frame(Forecasted = rep(5.2, 8), 
                            Hindcasted = rep(1.5, 8),
                            shelf = c(rep("West", 4), rep("East" , 4)),
                            monthf = c(rep(unique(ploto$monthf), 2)),
                            label = c(paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 6], ploto$Hindcasted[ploto$shelf == "West" & ploto$month == 6]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 7], ploto$Hindcasted[ploto$shelf == "West" & ploto$month == 7]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 8], ploto$Hindcasted[ploto$shelf == "West" & ploto$month == 8]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 9], ploto$Hindcasted[ploto$shelf == "West" & ploto$month == 9]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 6], ploto$Hindcasted[ploto$shelf == "East" & ploto$month == 6]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 7], ploto$Hindcasted[ploto$shelf == "East" & ploto$month == 7]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 8], ploto$Hindcasted[ploto$shelf == "East" & ploto$month == 8]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 9], ploto$Hindcasted[ploto$shelf == "East" & ploto$month == 9]), 2), sep = "")
                            ))
ann_text_pfsm$monthf <- month.abb[ann_text_pfsm$monthf]
ann_text_pfsm$monthf <- parse_factor(ann_text_pfsm$monthf, levels = month.abb)
ann_text_pfsm$shelf <- factor(ann_text_pfsm$shelf, levels = c("West", "East"))
(round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 6], ploto$Hindcasted[ploto$shelf == "East" & ploto$month == 6]), 2))
(round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 7], ploto$Hindcasted[ploto$shelf == "East" & ploto$month == 7]), 2))
(round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 8], ploto$Hindcasted[ploto$shelf == "East" & ploto$month == 8]), 2))
(round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 9], ploto$Hindcasted[ploto$shelf == "East" & ploto$month == 9]), 2))

(round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 6], ploto$Hindcasted[ploto$shelf == "West" & ploto$month == 6]), 2))
(round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 7], ploto$Hindcasted[ploto$shelf == "West" & ploto$month == 7]), 2))
(round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 8], ploto$Hindcasted[ploto$shelf == "West" & ploto$month == 8]), 2))
(round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 9], ploto$Hindcasted[ploto$shelf == "West" & ploto$month == 9]), 2))

hind_fore_plot_two_shelves_monthly <- ggplot(ploto, aes(y = Hindcasted, Forecasted))+
  geom_point(shape = 21, alpha = 1)+
  facet_grid(shelf ~ monthf)+
  geom_abline(slope = 1, color = "grey50", size = 0.3)+
  theme_bw()+
  theme(legend.position = "top")+
  scale_color_manual(values = Palette_months())+
  coord_fixed(ratio = 1, xlim = c(1, 6), ylim = c(1, 6), expand = TRUE)+
  labs(color = "", x = expression(Forecasted~BWDO~(mg/L)), y = expression(Hindcasted~BWDO~(mg/L)))+
  geom_text(data = ann_text_pfsm, label = ann_text_pfsm$label, parse = T, color = "black", size = 3)+
  theme(strip.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,-10,-10))
hind_fore_plot_two_shelves_monthly
#####################################################################################
## June Forecasted - Observed
ann_text_hosm <- data.frame(Observed= rep(1.5, 8),
                            Forecasted = rep(4.9, 8), 
                            shelf = c(rep("West", 4), rep("East" , 4)),
                            monthf = c(rep(unique(ploto$monthf), 2)),
                            label = c(paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 6], ploto$Observed[ploto$shelf == "West" & ploto$month == 6]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 7], ploto$Observed[ploto$shelf == "West" & ploto$month == 7]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 8], ploto$Observed[ploto$shelf == "West" & ploto$month == 8]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "West" & ploto$month == 9], ploto$Observed[ploto$shelf == "West" & ploto$month == 9]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 6], ploto$Observed[ploto$shelf == "East" & ploto$month == 6]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 7], ploto$Observed[ploto$shelf == "East" & ploto$month == 7]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 8], ploto$Observed[ploto$shelf == "East" & ploto$month == 8]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(ploto$Forecasted[ploto$shelf == "East" & ploto$month == 9], ploto$Observed[ploto$shelf == "East" & ploto$month == 9]), 2), sep = "")
                            ))
ann_text_hosm$monthf <- month.abb[ann_text_hosm$monthf]
ann_text_hosm$monthf <- parse_factor(ann_text_hosm$monthf, levels = month.abb)
ann_text_hosm$shelf <- factor(ann_text_hosm$shelf, levels = c("West", "East"))

fore_obs_plot_two_shelves_monthly <- ggplot(ploto, aes(x = Forecasted, Observed))+
  geom_point(shape = 21, alpha = 1)+
  facet_grid(shelf ~ monthf)+
  geom_abline(slope = 1, color = "grey50", size = 0.3)+
  theme_bw()+
  theme(legend.position = "top")+
  scale_color_manual(values = Palette_months())+
  coord_fixed(ratio = 1, xlim = c(1, 6), ylim = c(1, 6), expand = TRUE)+
  labs(color = "", y = expression(Observed~BWDO~(mg/L)), x = expression(Forecasted~BWDO~(mg/L)))+
  geom_text(data = ann_text_hosm, label = ann_text_hosm$label, parse = T, color = "black", size = 3)+
  theme(strip.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -10, -10))
fore_obs_plot_two_shelves_monthly
#####################################################################################
######################### This one
glimpse(ploto)
ploto_m <- ploto %>% dplyr::select(-month, -year, -monthf) %>%
  melt(., id.vars = c("Date", "shelf", "Forecasted"))
glimpse(ploto_m)
labels_fho <- c(expression(Hindcasted~(mg/L)),
                expression(Observed~(mg/L)))
ploto_m$variable <- factor(ploto_m$variable, 
                                          levels = c("Hindcasted", "Observed"), 
                                          ordered = TRUE, labels = labels_fho)
ploto_m$shelf <- factor(ploto_m$shelf, 
                                       levels = c("West", "East"), 
                                       ordered = TRUE)
#Rsquared values
(ofe <- format(round(E(ploto$Forecasted[ploto$shelf == "East"], ploto$Observed[ploto$shelf == "East"]), 2), nsmall = 2))
(ofw <- format(round(E(ploto$Forecasted[ploto$shelf == "West"], ploto$Observed[ploto$shelf == "West"]), 2), nsmall = 2))
(hfe <- format(round(E(ploto$Forecasted[ploto$shelf == "East"], ploto$Hindcasted[ploto$shelf == "East"]), 2), nsmall = 2))
(hfw <- format(round(E(ploto$Forecasted[ploto$shelf == "West"], ploto$Hindcasted[ploto$shelf == "West"]), 2), nsmall = 2))

ann_text_r2_fho_shelves <- data.frame(value = rep(1.7, 4),
                                      Forecasted = rep(5.1, 4),
                                      variable = c("Hindcasted ~ (mg/L)",
                                                   "Hindcasted ~ (mg/L)",
                                                   "Observed ~ (mg/L)",
                                                   "Observed ~ (mg/L)"),
                                      shelf = c("West", "East", "West", "East"),
                                      #label = c(paste("italic(R) ^ 2 == ", hfw, sep = ""),
                                      #          paste("italic(R) ^ 2 == ", hfe, sep = ""),
                                      #          paste("italic(R) ^ 2 == ", ofw, sep = ""),
                                      #          paste("italic(R) ^ 2 == ", ofe, sep = "")))
                                      label = c(hfw, hfe, ofw, ofe))

                                      
ann_text_r2_fho_shelves$shelf <- factor(ann_text_r2_fho_shelves$shelf, 
                                        levels = c("West", "East"), 
                                        ordered = TRUE)
ploto_m <- mutate(ploto_m, month = month(Date))
ploto_m$monthf <- month.abb[ploto_m$month]
ploto_m$monthf <- parse_factor((ploto_m$monthf), levels = month.abb)
hfo_bwdo <- ggplot(ploto_m, aes(y = value, x = Forecasted, color = monthf))+ #, color = monthf
  geom_point(shape = 21)+
  facet_rep_grid(variable ~ shelf, switch = "y", labeller = label_parsed, repeat.tick.labels = F)+
  geom_abline(slope = 1)+
  theme_few()+
  #geom_text(data = ann_text_r2_fho_shelves, 
            #label = sprintf(ann_text_r2_fho_shelves$label), 
  #          parse = T, color = "black", size = 3)+
  geom_text(data = ann_text_r2_fho_shelves[ann_text_r2_fho_shelves$shelf == "West" & ann_text_r2_fho_shelves$variable == "Hindcasted ~ (mg/L)",], 
            label = paste0("italic(R) ^ 2 == ", deparse(ann_text_r2_fho_shelves$label[1])), parse = T, color = "black", size = 3)+
  geom_text(data = ann_text_r2_fho_shelves[ann_text_r2_fho_shelves$shelf == "East" & ann_text_r2_fho_shelves$variable == "Hindcasted ~ (mg/L)",], 
            label = paste0("italic(R) ^ 2 == ", deparse(ann_text_r2_fho_shelves$label[2])), parse = T, color = "black", size = 3)+
  geom_text(data = ann_text_r2_fho_shelves[ann_text_r2_fho_shelves$shelf == "West" & ann_text_r2_fho_shelves$variable == "Observed ~ (mg/L)",], 
            label = paste0("italic(R) ^ 2 == ", deparse(ann_text_r2_fho_shelves$label[3])), parse = T, color = "black", size = 3)+
  geom_text(data = ann_text_r2_fho_shelves[ann_text_r2_fho_shelves$shelf == "East" & ann_text_r2_fho_shelves$variable == "Observed ~ (mg/L)",], 
            label = paste0("italic(R) ^ 2 == ", deparse(ann_text_r2_fho_shelves$label[4])), parse = T, color = "black", size = 3)+
  #scale_color_manual(values = Palette_months())+
  scale_color_brewer(palette = "Spectral")+
  labs(x = expression(Forecasted~BWDO~(mg/L)), y = "", color = "")+
  #scale_x_continuous(expand = c(0.05, 0.05), limits = c(0, 25000), breaks = c(0, 10000, 20000))+
  #scale_y_continuous(expand = c(0.05, 0.05), limits = c(0, 25000), breaks = c(0, 10000, 20000))+
  coord_fixed(ratio = 1, xlim = c(1.5, 6), ylim = c(1.5, 6), expand = TRUE)+
  theme(legend.position = "top",
        legend.text = element_text(size = 8),
        strip.placement = "outside",
        strip.text = element_text(size = 10), # text size inside strip (y-axis title)
        axis.title.x = element_text(size = 10), # text size inside strip (x-axis title) 
        axis.text = element_text(size = 8), # text size axis values
        legend.margin = margin(0,0,0,0),
        panel.spacing = unit(0.5, "lines"),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"), #top,right,bot, left
        #axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.9),
        legend.box.margin = margin(-7,15,-5,-7))
hfo_bwdo
date_now <- Sys.Date()
png(filename = paste("./graphics for manuscript/bwdo_hind_obs_vs_fore_two_shelves_", date_now, ".png", sep = ''),
    units = "cm", 
    width = 14, 
    height = 12, 
    #pointsize=12, 
    res = 500)
print(hfo_bwdo)
dev.off()
