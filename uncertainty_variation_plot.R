if (!require("pacman")) install.packages("pacman"); 
pacman::p_load(magrittr, tidyverse, lubridate, MHadaptive, IDPmisc, 
               reshape2, ggfortify, tictoc, ggthemes, faraway, rstudioapi, ggpubr) 
select <- dplyr::select
#### Get directory - location of the current R file and set it as basepath
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))
getwd()
#####################################################################################
# All functions here
source("all_functions.R")
#####################################################################################
# Load forecast, hindcast simulations
out_plot <- data.frame()
path <- paste0(getwd(), "/data/all_sim_09-09_paral/")
years_of_interest <- c(1985:2016)
for (i in seq_along(years_of_interest)) {
  yy <- years_of_interest[i]
  temp_data <- readRDS(file = paste(path, "ha_old_", yy, ".rds", sep = ""))
  out_plot <- bind_rows(out_plot, temp_data)
}
dplyr::glimpse(out_plot)
ci_out <- out_plot %>% 
  dplyr::select(Date, ha_ParHM_ba_mean, `ha_ParHM_ba_97.5%`, `ha_ParHM_ba_2.5%`, ha_ParHM_ba_mean)
ci_out$monthday <- format(as.Date(ci_out$Date), "%m-%d")
ci_out <- ci_out %>%  
  mutate(ha_97.5 = `ha_ParHM_ba_97.5%` / ha_ParHM_ba_mean,
         ha_2.5 = `ha_ParHM_ba_2.5%` / ha_ParHM_ba_mean)
glimpse(ci_out)
#####################################################################################

############## 
combined_plot_5 <- bind_cols(ci_out %>% 
                               mutate(m = ha_ParHM_ba_mean,
                                      iqr = `ha_ParHM_ba_97.5%` - `ha_ParHM_ba_2.5%`) %>%
                               dplyr::select(monthday, m, iqr) %>%
                               melt(., id.vars = "monthday"),
                             ci_out %>%  
                               mutate(norm_iqr = ha_97.5 - ha_2.5) %>%
                               dplyr::select(norm_iqr) %>%
                               bind_rows(., tibble(norm_iqr = rep(NA, nrow(ci_out)))))
glimpse(combined_plot_5)
m <- combined_plot_5 %>% 
  group_by(monthday, variable) %>% 
  summarise(value = mean(value, na.rm = T))
glimpse(m)
ratio <- max(na.omit(combined_plot_5$norm_iqr)) / max(na.omit(m$value[m$variable == "iqr"]))
iqr_all_plot <- ggplot() +
  geom_boxplot(data = combined_plot_5, outlier.shape = NA, aes(x = monthday, y = norm_iqr), lwd = 0.3, fatten = 1)+ 
  geom_line(data = m, aes(x = monthday, y = value * ratio , group = variable, linetype = variable), color = "red", size = 0.9)+
  theme_few()+
  scale_y_continuous(#labels = ,
    sec.axis = sec_axis(~./ratio/ 10 ^ 4, name = expression(HA~(10^{4}~km^{2})))
  )+
  labs(y = "Normalized 95% HA IQR", x = "", linetype = "")+
  scale_x_discrete(expand = c(0, 0), 
                   breaks = c("06-01", "07-01", "08-01", "09-01"),
                   labels = c("Jun", 'Jul', 'Aug', 'Sep'))+
  #scale_y_continuous(limits = c(-60, 60))+
  scale_linetype_manual(values = c("dashed","solid"), labels = c(m = "Mean HA", iqr = "Mean 95% HA IQR"))+
  annotate("text", label = "A", size = 4, x = -Inf, y = Inf, vjust=1.2, hjust=-0.3)+
  theme(legend.position = "top",
        legend.key.width = unit(1,"cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.text = element_text(color = "red", size = 8),
        axis.title = element_text(size = 8),  
        axis.title.y.right = element_text(color = "red"),
        axis.text.y.right = element_text(color = "red"),
        axis.line.y.right = element_line(color = "red"),
        axis.ticks.y.right = element_line(color = "red"),
        axis.text = element_text(size = 8), # text size axis values
        strip.placement = "outside",
        strip.text = element_text(size = 8, vjust=-0.5), # text size inside strip
        plot.margin = unit(c(0.5, 0.1, 0.2, 0.2),"cm")) #top,right,bot, left
iqr_all_plot
### Winds plot
InpDai <- read_csv(file="./data/InpDaiGomMod_mar18v2.csv")
InpDai <- InpDai[, 2:ncol(InpDai)]
InpDai <- InpDai %>% 
  mutate(year = year(outp_tim),
         month = month(outp_tim), 
         monthday = format(as.Date(outp_tim), "%m-%d")) %>% 
  dplyr::rename(Date = outp_tim)
str(InpDai)
winds_plot_data <- ggplot(InpDai, aes(monthday, smr14wi2W))+
  geom_boxplot(fill = "grey80", outlier.shape = NA, lwd = 0.3, fatten = 1)+
  theme_few()+
  labs(y = expression(Wind~speed^{2}~(m/s)^{2}), x = "")+
  scale_x_discrete(expand = c(0, 0), 
                   breaks = c("06-01", "07-01", "08-01", "09-01"),
                   labels = c("Jun", 'Jul', 'Aug', 'Sep'))+
  #scale_y_continuous(limits = c(-60, 60))+
  annotate("text", label = "B", size = 4, x = -Inf, y = Inf, vjust=1.2, hjust=-0.3)+
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.y.right = element_text(color = "red"),
        axis.text = element_text(size = 8), # text size axis values
        strip.placement = "outside",
        axis.line.y.right = element_line(color = "red"),
        strip.text = element_text(size = 8, vjust=-0.5), # text size inside strip
        plot.margin = unit(c(-0.4, 0.1, 0.2, 0.2),"cm")) #top,right,bot, left
winds_plot_data
ggarrange(iqr_all_plot, winds_plot_data, ncol = 1, align = "v", heights = c(1.5, 1))
(sysdate <- Sys.Date())
png(filename = paste("./graphics for manuscript/ci_real_mean_norm_ws_", sysdate, '.png', sep = ''), 
    units = "cm", 
    width = 16, 
    height = 9, 
    #pointsize=12, 
    res = 500)
print(ggarrange(iqr_all_plot, winds_plot_data, ncol = 1, align = "v", heights = c(1.5, 1)))
dev.off()
