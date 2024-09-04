if (!require("pacman")) install.packages("pacman"); 
pacman::p_load(adaptMCMC, TTR, truncnorm, tidyverse, lubridate, suncalc, data.table, deSolve,
               zoo, faraway, gridExtra, reshape2, DataCombine, leaps, gdata, ggthemes,
               cowplot, scales, GGally, odin, gtable, rstudioapi, matrixStats, hydroGOF, EGRET,
               rnoaa, tsibble, feasts, rsample, glmnet, ggpubr, sjPlot)
select <- dplyr::select
#### Get directory - location of the current R file and set it as basepath
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))
getwd()
#####################################################################################
# All functions here
source("all_functions.R")
#####################################################################################
## Read Monthly
dailyseq <- seq(as.Date("1967/10/1"), as.Date("2016/9/1"), "day")
monthlyseq <- data.frame(Date = dailyseq) %>% 
  mutate(year = year(Date), month = month(Date)) %>%
  group_by(year, month) %>%
  summarise(daysinmonth = as.numeric(n())) # also counting number of days in a month
gommonthly <- readxl::read_xlsx("./data/Gulf-Monthly-2016.xlsx", skip = 40, col_names = F)
gomm <- bind_cols(monthlyseq, gommonthly %>% select(AtQm = 6, AtNOx = 57, AtTKN = 61, AtNH3 = 65, 
                                                    MsQm = 3, MsNOx = 32, MsTKN = 37, MsNH3 =  41))
gom <- gomm %>% 
  mutate(AtNbio = (AtNOx + AtNH3 + (AtTKN - AtNH3) * .12), # tonns / second
         MsNbio = (MsNOx + MsNH3 + (MsTKN - MsNH3) * .12),
         AtCm = AtNbio / AtQm,
         MsCm = MsNbio / MsQm) %>% # this is in mega gramms (gamms * 10^-6) because loading was in metric tonns
  select(year, month, AtQ = AtQm, AtN = AtNbio, MsQ = MsQm, MsN = MsNbio) %>%
  filter(year > 1979)
ggplot(gom, aes(AtN/AtQ, MsN/MsQ))+geom_point()
head(gom)
glimpse(gom)
RivDaim <- gom
cor(gom$AtQ, gom$AtN);
cor(gom$MsQ, gom$MsN)


Qy <- gom %>% 
  group_by(year) %>%
  filter(month %in% c(6,7,8,9)) %>% 
  summarise(Qmedian = median(MsN, na.rm = T))

Qstat <- quantile(Qy$Qmedian, probs = c(0.333, 0.666, 1))
Qy$hytype <- NA
for (i in 1:nrow(Qy)) {
  if (Qy$Qmedian[i] <= Qstat[1]) {
    Qy$hytype[i] <- "Dry"
  } else if (Qy$Qmedian[i] <= Qstat[2]) {
    Qy$hytype[i] <- "Normal" 
  } else {
    Qy$hytype[i] <- "Wet"
  }
}
Qy %>% filter(year > 1984) %>% group_by(hytype) %>% summarize(n())
## 
##
## Read precipitation and temperature
precip <- read_csv("./data/prism_prec_monthly.csv") %>% 
  mutate(p = p * 2.54) # convert to cm
glimpse(precip)
temper <- read_csv("./data/prism_temp_monthly.csv")
predictors <- left_join(precip, temper, by = c("year", "month")) #%>% select(-am_P, -jfm_T)
win <- predictors %>% 
  filter(month %in% c(1:4)) %>%
  group_by(year) %>%
  summarise(P_1234 = mean(p, na.rm = T),
            T_1234 = mean(t, na.rm = T))  
may <- predictors %>% 
  filter(month %in% c(5)) %>%
  group_by(year) %>%
  summarise(P_5 = mean(p, na.rm = T),
            T_5 = mean(t, na.rm = T))
predictorss <- left_join(win, may, by = "year")
##
cor(RivDaim %>% filter(month == 5) %>% pull(AtQ), predictorss$P_5)
cor(RivDaim %>% filter(month == 5) %>% pull(MsQ), predictorss$P_5)

cor(RivDaim %>% filter(month == 7) %>% pull(AtQ), predictorss$P_5)
cor(RivDaim %>% filter(month == 7) %>% pull(MsQ), predictorss$P_5)
cor(RivDaim %>% filter(month == 6) %>% pull(AtQ), predictorss$P_5)
cor(RivDaim %>% filter(month == 6) %>% pull(MsQ), predictorss$P_5)
cor(RivDaim %>% filter(month == 8) %>% pull(AtQ), predictorss$P_5)
cor(RivDaim %>% filter(month == 8) %>% pull(MsQ), predictorss$P_5)


months_to_predict <- c(6, 7, 8, 9)
respnames <- c("AtQ", "AtN", "MsQ", "MsN")
mod_list <- list()
for (i in seq_along(respnames)) {
  varname <- respnames[i]
  del <- RivDaim %>%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    filter(month %in% c(1:9)) %>%
    select(year, month, varname) %>% 
    dcast(., year ~ month, value.var = varname) %>% 
    mutate(., `1234` = (`1` + `2` + `3` + `4`) / 4) %>% 
    select(-c(`1`, `2`, `3`, `4`))
    colnames(del)[2:ncol(del)] <- paste(varname, colnames(del)[2:ncol(del)], sep = "_")     
    assign(paste(varname, sep = ""), del)
  savelist <- list()
  for (j in seq_along(months_to_predict)) {
  var_to_model <- paste(varname, months_to_predict[j], sep = "_")
  alldat <- left_join(del, predictorss, by = "year") %>% 
    select(response = var_to_model, colnames(del)[c(2, 7)], 
           colnames(predictorss[2:ncol(predictorss)]))
    rs <- regsubsets(sqrt(response) ~ ., alldat, nvmax = 2)
    summary(rs)
    plot(rs)
    plot(rs, scale = "r2")
    bic <- summary(rs)$bic
    bic.idx <- which(bic == min(bic))
    sel <- summary(rs)$which[bic.idx, ]
    names_sel <- names(sel[sel == T])[-1]
    preds <- colnames(alldat %>% dplyr::select(names_sel))
    mdl <- lm(sqrt(response) ~ ., alldat %>% dplyr::select(response, names_sel))
    summary(mdl)
    vif(mdl)
    #plot(predict(mdl), sqrt(mdl$model$response), xlab = "Predicted", ylab = "Observed"); abline(a = 0, b = 1)
    savelist[[j]] <- list(response = var_to_model, model = mdl, vif = vif(mdl), 
                          rsq = summary(mdl)$adj.r.squared, preds)
    }
  mod_list <- append(mod_list, savelist)
}

summary(mod_list)
#save(mod_list, file = "./data/AMQL_models_sqrt_prism_05-11.RData")
#load("./data/AMQL_models_sqrt_prism_05-11.RData")
## Some analysis for weekly meetings
pr <- data.frame(preds = matrix(NA, nrow = 16, ncol = 1))
for (i in seq_along(mod_list)) {
  pr[i, ] <- paste0(unlist(mod_list[[i]][[5]]), collapse = ", ")
} 
tibble(response = sapply(mod_list, function(x) x[[1]]),
       predictors = pr$preds,
       adj.rsquared = round(sapply(mod_list, function(x) x[[4]]), 2)
       )

sapply(mod_list, function(x) x[[3]])
sapply(mod_list, function(x) x[[1]])
for (i in seq_along(mod_list)) {
  print(paste(mod_list[[i]]$response))
  (sumary(mod_list[[i]]$model))
}
mod_list[[1]]$model$coefficients
sapply(mod_list, function(x) (x[[5]]))
for (i in seq_along(mod_list)) { 
  print(paste(mod_list[[i]]$response))
  print((mod_list[[i]]$model$coefficients))
  }
##
## List of coefficients
(newoldcoefs <- bind_cols(response = sapply(mod_list, function(x) x[[1]]), 
                        new_coefs = rep(NA, length(mod_list)),
                        old_coefs = rep(NA, length(mod_list))))
for (i in seq_along(mod_list)) {
 if (any(names(coef(mod_list[[i]]$model)) == "P_5")) {
   
   newoldcoefs$new_coefs[i] <- round((coef(mod_list[[i]]$model)[names(coef(mod_list[[i]]$model)) 
                                                                == "P_5"]), 2)
   newoldcoefs$old_coefs[i] <- newoldcoefs$new_coefs[i] * 2.54
   }  
}
newoldcoefs

length(mod_list[[1]]$model$coefficients)

## Check for auto correlation. plotting residuals
for (i in seq_along(mod_list)) {
  n <- length(residuals(mod_list[[i]][[2]]))
  sumary(lm(tail(residuals(mod_list[[i]][[2]]), n-1) ~ head(residuals(mod_list[[i]][[2]]), n-1) -1))
}
for (i in seq_along(mod_list)) {
  n <- length(residuals(mod_list[[i]][[2]]))
  print(cor(tail(residuals(mod_list[[i]][[2]]), n-1),head(residuals(mod_list[[i]][[2]]), n-1)))
}

# That is the axis and plot names for images
units_exp_observed <- c(rep("Observed~(sqrt(m^{3}/s))", 4), rep("Observed~(sqrt(T/mo))", 4), 
                        rep("Observed~(sqrt(m^{3}/s))", 4), rep("Observed~(sqrt(T/mo))", 4))
units_exp_predicted <- c(rep("Predicted~(sqrt(m^{3}/s))", 4), rep("Predicted~(sqrt(T/mo))", 4), 
                         rep("Predicted~(sqrt(m^{3}/s))", 4), rep("Predicted~(sqrt(T/mo))", 4))
response_name <- c("sqrt(italic(Q)[italic(A)][6])", "sqrt(italic(Q)[italic(A)][7])",
                   "sqrt(italic(Q)[italic(A)][8])", "sqrt(italic(Q)[italic(A)][9])",
                   "sqrt(italic(L)[italic(A)][6])", "sqrt(italic(L)[italic(A)][7])",
                   "sqrt(italic(L)[italic(A)][8])", "sqrt(italic(L)[italic(A)][9])",
                   "sqrt(italic(Q)[italic(M)][6])", "sqrt(italic(Q)[italic(M)][7])",
                   "sqrt(italic(Q)[italic(M)][8])", "sqrt(italic(Q)[italic(M)][9])",
                   "sqrt(italic(L)[italic(M)][6])", "sqrt(italic(L)[italic(M)][7])",
                   "sqrt(italic(L)[italic(M)][8])", "sqrt(italic(L)[italic(M)][9])"
)

adjusted_rsq <- round(sapply(mod_list, function(x) x[[4]]), 2)

expressions_for_plots <- cbind.data.frame(name = response_name, units_obs = units_exp_observed,
                                          units_pred = units_exp_predicted)

plot_list <- list()
time_now <- Sys.Date()
pdf(file = paste("./data/regressions_sqrt_", time_now, ".pdf", sep=''), onefile = T)
## This creates an empty table
table_df_list <- list()
for (i in 1:length(mod_list)) {
  tt <- (mod_list[[i]][[2]])
  #par(mfrow = c(2, 2))
  ## This is for table
  temp_table <- cbind.data.frame(predicted = (tt$fitted.values) ^ 2,
                                 observed = (tt$model$`sqrt(response)`) ^ 2)
  colnames(temp_table) <- c(paste0('Predicted ', mod_list[[i]][[1]]), paste0('Observed ', mod_list[[i]][[1]]))
  table_df_list[[i]] <- temp_table
  ## This is for plots
  ttt <- cbind.data.frame(predicted = tt$fitted.values,
                          observed = tt$model$`sqrt(response)`)
  name <- expressions_for_plots[i, 1]
  units_obs <- expressions_for_plots[i, 2]
  units_pred <- expressions_for_plots[i, 3]
  max_scale <- max(max(ttt$observed), max(ttt$predicted))
  min_scale <- min(min(ttt$observed), min(ttt$predicted))
  
  plotit <- ggplot(ttt, aes(predicted, observed))+
    geom_point(shape = 1)+
    theme_few()+
    geom_abline(slope = 1)+
    ggtitle(parse(text = name))+
    labs(x = parse(text = units_pred),
         y = parse(text = units_obs))+
    scale_x_continuous(limits = c(min_scale-5, max_scale+5))+
    scale_y_continuous(limits = c(min_scale-5, max_scale+5))+
    theme(plot.title = element_text(size = 10, hjust = 0.5),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          aspect.ratio = 1)
  plotit
  print(plotit)
  plot_list[[i]] <- plotit
}
dev.off()

table_df <- do.call(cbind, table_df_list)
table_df <- bind_cols(Year = 1980:2016, table_df)
## Save table
write_csv(table_df, "./graphics for manuscript/monthly_riverine.csv")

for (i in c(1,5,9,13)) {
  plott <- ggarrange(plot_list[[i]], plot_list[[(i+1)]], 
                     plot_list[[(i+2)]], plot_list[[(i+3)]], 
                     ncol = 2, nrow = 2)
  ggsave(
    filename = paste("./graphics for manuscript/regressions_", i, "_", time_now, ".png", sep = ""),
    plot = plott,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 15,
    height = 15,
    units = c("cm"),
    dpi = 500)
}

  