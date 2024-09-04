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
# load updated inputs
InpDai <- read_csv(file="./data/InpDaiGomMod_mar18v2.csv")
InpDai <- InpDai[, 2:ncol(InpDai)]
InpDai <- InpDai %>% 
  mutate(year = year(outp_tim),
         month = month(outp_tim)) %>% 
  dplyr::rename(Date = outp_tim)
str(InpDai)
days_of_interest <- InpDai %>% dplyr::select(Date)
#
# load output data (from cruises and geostat model)
Out.dat <- read_csv(file = "./data/GoM_cruise_sep18.csv", skip = 1)
Out.dat$Date <- as.Date(Out.dat$Date, format = "%m/%d/%Y")
str(Out.dat)
Out.dat.all <- left_join(days_of_interest, Out.dat, by = "Date") %>% 
  dplyr::select(Date, obsOw = Mean_4,   obsOe = Mean_3, # define observations 
                errOw = STDEV_4,  errOe = STDEV_3) # define sdv geostat errors 
##  
obs_index <- which(is.na(Out.dat.all$obsOe) == F) # vector with dates observed
# ---------------------------------------------------------------------------------------------
# add Os observations (actually load surface wT and S data and transform them into Os) --------
# ---------------------------------------------------------------------------------------------
## not changed
wT0 <- read_csv(file="./data/shelf_salinity_temp_apr18.csv", na = "#N/A")
glimpse(wT0)
wT <- wT0[,c("Date","Surface Temp","Surface Salinity","Surface Temp_1","Surface Salinity_1")]
wT$Date <- as.Date(wT$Date, format = "%m/%d/%Y")
glimpse(wT)
wT <- wT %>% arrange(Date) # orders rows chronologically
wT <- as.data.frame(wT) # otherwise the tricky code doesn't work =\
Out.dat <- as.data.frame(Out.dat)
# replace these dates with exactly the same dates as the cruises (and take care of duplicates)
for (i in 1:nrow(wT)){
  wT[i,1] <- Out.dat[which(abs(wT[i,1]-Out.dat[,1]) == min(abs(wT[i,1] - Out.dat[,1]))),1][1]
}
# estimate wT and S, then Os (in the surface)
Tsrf.dy <- data.frame(time = InpDai[obs_index, 1], Tsrf = NA)
for (i in 1:nrow(InpDai[obs_index, ]))
{
  ndx.rw <- which(wT[, 1] == Tsrf.dy[i, 1]) 
  if(length(ndx.rw)>0) 
  {
    Tsrf.dy[i, 2] = mean(as.matrix(wT[ndx.rw,c("Surface Temp","Surface Temp_1")]), na.rm = T)
  }
}
Ssrf.dy <- data.frame(time = InpDai[obs_index, 1], Tsrf = NA)
for (i in 1:nrow(InpDai[obs_index, ])) {
  ndx.rw <- which(wT[,1] == Ssrf.dy[i,1]) 
  if(length(ndx.rw)>0) {
    Ssrf.dy[i,2] = mean(as.matrix(wT[ndx.rw,c("Surface Salinity","Surface Salinity_1")]),na.rm = T)
  }
}
Tsrf.dy[,2]=Tsrf.dy[,2]+ 273.15
o_ss = exp((-139.34411+1.575701*1e5/Tsrf.dy[,2]-6.642308*1e7/Tsrf.dy[,2]^2+1.243800*1e10/Tsrf.dy[,2]^3-8.621949*1e11/Tsrf.dy[,2]^4)-
             Ssrf.dy[,2]*(1.7674*1e-2 - 1.0754*1e1/Tsrf.dy[,2]+2.1407*1e3/Tsrf.dy[,2]^2)) # saturation concentration of dissolved oxygen in saltwater at 1 atm (mg L-1)
modOs <- data.frame(time = InpDai[obs_index,1], satO = o_ss)
#
# use month-dependent constant Os
# ----------------------------------------------------------
meOs <- cbind(6:9, NA)
for (i in 6:9)
{ ndx <- which(month(modOs[,1])==i)
modOs[ndx,2] <- mean(modOs[ndx,2],na.rm = T)
meOs[i-5,2]<- mean(modOs[ndx,2],na.rm = T)
}
colnames(meOs) <- c("month", "Os")
meOs <- as.data.frame(meOs)
# ---------------------------------------------------------------------------------------------
# calculate month-specific bwT
wT0 <- read_csv(file="./data/shelf_salinity_temp_apr18.csv", skip=0)
mon_bwT <- wT0 %>% dplyr::select(cru_date="Date",bwT_W="30m temp",bwT_E="30m temp_1")%>% 
  mutate(cru_date=mdy(cru_date),bwT_W=as.numeric(bwT_W),bwT_E=as.numeric(bwT_E)) %>% 
  transmute(cru_date,bwT_shlf=rowMeans(cbind(bwT_W,bwT_E), na.rm = T)) %>% 
  group_by(month(cru_date)) %>% summarize(month_avg_bwT=mean(bwT_shlf, na.rm = T)) %>% 
  filter(.[[1]]>5)
ref_bwT <- mon_bwT %>% pull(2) %>% mean
colnames(mon_bwT) <- c("month", "bWt")
# ------------------------------------------------------------
# Calculate winter loads
# ------------------------------------------------------------
RivDai <- read.csv(file = "./data/daiAtMs_feb18.csv", header = T, skip = 0, sep = ",")
RivDai[, 1] <- seq(as.Date("1980/01/1"), as.Date("2016/09/30"), "days")
RivDai[, c(3,5)] <- RivDai[ , c(3, 5)] * 30.5 * 86400 # conversion to /mo
colnames(RivDai)[c(3, 5)] <-c("At_bioN_Mgpmo", "Ms_bioN_Mgpmo")
OldLoaDai <- data.frame(Date = days_of_interest) %>% add_column(LoNoMa = NA)

for (i in 1:length(OldLoaDai$Date)) 
{
  prim_nov <- as.Date(paste("November 1,", year(OldLoaDai[i, 1]) - 1, sep = " "), format = '%B %d, %Y')
  ulti_mar <- as.Date(paste("March 31,",   year(OldLoaDai[i, 1]),     sep = " "), format = '%B %d, %Y')
  perNovMar <- seq(prim_nov, ulti_mar, "days")
  OldLoaDai$LoNoMa[i] <- RivDai %>% as_tibble %>% filter(time %in% perNovMar) %>%  
    transmute(date_day=time,tot_loa=rowSums(cbind(At_bioN_Mgpmo,Ms_bioN_Mgpmo))) %>% 
    pull(tot_loa) %>% mean # average total monthly load (now Atc&Miss are summed)
}
avgNwint <- mean(unique(OldLoaDai$LoNoMa))
avgNwint
###
glimpse(InpDai)
InpDaif <- left_join(InpDai, OldLoaDai, by = "Date") %>%
  left_join(., meOs, by = "month") %>% 
  left_join(., mon_bwT, by = "month")
const_input <- left_join(mon_bwT, meOs, by = "month") 

# get the variables of interest
vars_of_interest <- colnames(InpDaif %>% dplyr::select(-Date, -year, -month))
glimpse(InpDaif)
get_input_matrix <- function(data, variables) {
  for (i in 1:length(variables)) {
    var <- variables[i]
    data$monthday <- factor(format(data$Date, format = "%m-%d"))
    tt <- data %>% 
      dplyr::select(year, monthday, var) %>% 
      reshape2::dcast(., monthday ~ year) %>% 
      dplyr::select(-monthday)
    assign(paste(var), tt, envir = parent.frame())
  }
}
get_input_matrix(data = InpDaif, variables = vars_of_interest)
###
years_of_calibration <- c(1985:2016)
number_of_days <- 1:nrow(spgAtQ)    # June-Sep
###
par_oct18 <- readRDS(file = "./data/pars.rds")
par <- colMeans(par_oct18)
# fixed parameters, hardcoded in the model 
R_CN = 5.68; R_OC = 3.5; # ratio of carbon/nitr, oxyg/carbon
A.W = 48450 / 1000  #shelf area 1000km2, G(m2)
A.E = 14000 / 1000  #shelf area 1000km2, G(m2)
CS_O = 3 #reference DO for sediment oxygen demand
FML = 0.2 #fraction of MS River lost under all wind conditions
DOadj = 1 #adjustment to model lower layer DO, not BWDO
rQg = 3.2 # Ocean dilution factor
meanMsQ = mean(colMeans(spgMsQ))
## Create empty mat
bwdoE <- matrix(NA, nrow = length(number_of_days), ncol = length(years_of_calibration)); bwdoW <- bwdoE;
for (j in 1:length(years_of_calibration)) {
  for(i in seq_along(number_of_days)) {
    # East/Mississ bottom segm
    pngFe    <- 0.5 + par["bE"] * spgWwi[i, j] #eq 5, spring / WCOD related 
    pngFe    <- min(max(1e-2, pngFe), 1) # Fe constrained
    smrFe    <- 0.5 + par["bE"] * smrWwi[i, j] #eq 5, summer / reaeration related
    smrFe    <- min(max(1e-2, smrFe), 1)
    QsE      <- (smrAtQ[i, j] * 86400 / 10 ^ 9) * smrFe # for eq 6, summer freshwater input [Gm3/d]
    tauE     <- (smr14wi2W[i, j] / 360 + smr14wi2E[i, j] / 110) / (1 / 360 + 1 / 110) / 100
    kaE      <- par["b0Ka"] + par["b1Ka"] * tauE / QsE * A.E / 100 #eq. 6, summer reaerat [m/d]
    NinE     <- spgMsN[i, j] / 30.5 / 1000 * (1 - FML) * (1 - pngFe) + spgAtN[i, j] / 30.5 / 1000 * pngFe # load entering in spring (QrCrN+QuCuN) [Gg/d]
    QinE     <- (rQg * mean(colMeans(spgMsQ)) * 86400 / 10 ^ 9 + 
                   spgMsQ[i, j] * 86400 / 10 ^ 9 * (1 - FML)) * (1 - pngFe) + 
      spgAtQ[i, j] * 86400 / 10 ^ 9 * pngFe # [Gm3/d]
    CfxE  <- NinE * R_CN * par["vs"] / (A.E * par["vs"] + QinE) # gC/(m2*d)
    bwdoE[i, j] <- (kaE * Os[i, j] - CfxE * R_OC * par["R_PR"])/(kaE + (LoNoMa[i, j] / avgNwint) ^ (.5) *
                                                                   par["SD"] * 1.07 ^ (bWt[i, j] - ref_bwT) / CS_O) - DOadj
    
    # West/Atchafa bottom segm
    QsW      <- (smrAtQ[i, j] * 86400 / 10 ^ 9) * (1 - smrFe) # for eq 6, [Gm3/d]
    tauW     <- (smr14wi2W[i, j] / 130 + smr14wi2E[i, j] / 350)/ (1 / 130 + 1 / 350) / 100
    kaW      <- par["b0Ka"] + par["b1Ka"] * tauW / QsW * A.W / 100 # eq. 6, summer reaerat [m/d]
    NoutE    <- NinE / (A.E * par["vs"] + QinE) * QinE # reduced by setting [Gg/d]
    NinW     <- (NoutE + spgAtN[i, j] / 30.5 / 1000) * (1 - pngFe) # available loads moved westward
    QinW     <- (QinE + spgAtQ[i, j] * 86400 / 10^9) * (1 - pngFe) + 
      rQg * mean(colMeans(spgMsQ)) * 86400 / 10 ^ 9 * pngFe # [Gm3/d]
    CfxW     <- NinW * R_CN * par["vs"] / (A.W * par["vs"] + QinW)
    bwdoW[i, j] <- (kaW * Os[i, j] - CfxW * R_OC * par["R_PR"]) / (kaW + (LoNoMa[i, j] / avgNwint) ^ (.5) * 
                                                                     par["SD"] * 1.07 ^ (bWt[i, j] - ref_bwT) / CS_O) - DOadj #
  }
}
summary(bwdoE)
doe <- reshape2::melt(bwdoE)[, 3][obs_index]
dow <- reshape2::melt(bwdoW)[, 3][obs_index]
do <- append(dow, doe)
#####################################################################################
# for one year 
one_year_run <- function(data, params) {
  data <- data %>% add_column(bwdoE = NA, bwdoW = NA, bwdoW_e = NA, bwdoE_e = NA)
  for(i in seq_along(number_of_days)) {
    # East/Mississ bottom segm
    pngFe    <- 0.5 + params["bE"] * data$spgWwi[i] #eq 5, spring / WCOD related 
    pngFe    <- min(max(1e-2, pngFe), 1) # Fe constrained
    smrFe    <- 0.5 + params["bE"] * data$smrWwi[i] #eq 5, summer / reaeration related
    smrFe    <- min(max(1e-2, smrFe), 1)
    QsE      <- (data$smrAtQ[i] * 86400 / 10 ^ 9) * smrFe # for eq 6, summer freshwater input [Gm3/d]
    tauE     <- (data$smr14wi2W[i] / 360 + data$smr14wi2E[i] / 110) / (1 / 360 + 1 / 110) / 100
    kaE      <- params["b0Ka"] + params["b1Ka"] * tauE / QsE * A.E / 100 #eq. 6, summer reaerat [m/d]
    NinE     <- data$spgMsN[i] / 30.5 / 1000 * (1 - FML) * (1 - pngFe) + data$spgAtN[i] / 30.5 / 1000 * pngFe # load entering in spring (QrCrN+QuCuN) [Gg/d]
    QinE     <- (rQg * meanMsQ * 86400 / 10 ^ 9 + 
                   data$spgMsQ[i] * 86400 / 10 ^ 9 * (1 - FML)) * (1 - pngFe) + 
      data$spgAtQ[i] * 86400 / 10 ^ 9 * pngFe # [Gm3/d]
    CfxE  <- NinE * R_CN * params["vs"] / (A.E * params["vs"] + QinE) # gC/(m2*d)
    data$bwdoE[i] <- (kaE * data$Os[i] - CfxE * R_OC * params["R_PR"]) / (kaE + (data$LoNoMa[i] / avgNwint) ^ (.5) *
                                                                            params["SD"] * 1.07 ^ (data$bWt[i] - ref_bwT) / CS_O) - DOadj
    
    # West/Atchafa bottom segm
    QsW      <- (data$smrAtQ[i] * 86400 / 10 ^ 9) * (1 - smrFe) # for eq 6, [Gm3/d]
    tauW     <- (data$smr14wi2W[i] / 130 + data$smr14wi2E[i] / 350)/ (1 / 130 + 1 / 350) / 100
    kaW      <- params["b0Ka"] + params["b1Ka"] * tauW / QsW * A.W / 100 # eq. 6, summer reaerat [m/d]
    NoutE    <- NinE / (A.E * params["vs"] + QinE) * QinE # reduced by setting [Gg/d]
    NinW     <- (NoutE + data$spgAtN[i] / 30.5 / 1000) * (1 - pngFe) # available loads moved westward
    QinW     <- (QinE + data$spgAtQ[i] * 86400 / 10^9) * (1 - pngFe) + 
      rQg * meanMsQ * 86400 / 10 ^ 9 * pngFe # [Gm3/d]
    CfxW     <- NinW * R_CN * params["vs"] / (A.W * params["vs"] + QinW)
    data$bwdoW[i] <- (kaW * data$Os[i] - CfxW * R_OC * params["R_PR"]) / (kaW + (data$LoNoMa[i] / avgNwint) ^ (.5) * 
                                                                            params["SD"] * 1.07 ^ (data$bWt[i] - ref_bwT) / CS_O) - DOadj #
  }
  data$bwdoW_e <- data$bwdoW + rnorm(nrow(data), 0, params["sy.W"])
  data$bwdoE_e <- data$bwdoE + rnorm(nrow(data), 0, params["sy.E"])
  data <- data %>% dplyr::select(Date, bwdoW, bwdoW_e, bwdoE, bwdoE_e)
  return(data)
}

dd <- data.frame()
for (i in seq_along(years_of_calibration)) {
  tem <- one_year_run(data = InpDaif %>% filter(year == years_of_calibration[i]), params = par)
  dd <- bind_rows(dd, tem)
}
plot(tem$Date, tem$bwdoE)
de <- melt(bwdoE) %>% bind_cols(Date = InpDaif$Date, moe = .$value) %>% dplyr::select(Date, moe)
dw <- melt(bwdoW) %>% bind_cols(Date = InpDaif$Date, mow = .$value) %>% dplyr::select(Date, mow)
ded <- left_join(dd, de, by = "Date") %>% 
  left_join(., dw, by = "Date")
plot(ded$moe, ded$bwdoE); abline(0,1)
plot(ded$mow, ded$bwdoW); abline(0,1)
############ 
pred_bwdo <- dd
glimpse(pred_bwdo)
plot(pred_bwdo$bwdoW, pred_bwdo$bwdoW_e); abline(0,1)
plot(pred_bwdo$bwdoE, pred_bwdo$bwdoE_e); abline(0,1)
# get observed dissolved oxygen
glimpse(Out.dat.all)
pred_obs_bwdo <- left_join(pred_bwdo, Out.dat.all, by = "Date")
E(pred_obs_bwdo$bwdoW, pred_obs_bwdo$obsOw)
E(pred_obs_bwdo$bwdoE, pred_obs_bwdo$obsOe)
pred_obs_bwdo_fin <- bind_rows(pred_obs_bwdo %>% dplyr::select(Date, Hindcasted = bwdoE, 
                                           Observed = obsOe) %>%
                     mutate(shelf = "East"),
                     pred_obs_bwdo %>% dplyr::select(Date, Hindcasted = bwdoW, 
                                           Observed = obsOw) %>%
                     mutate(shelf = "West")) %>% 
  mutate(month = month(Date))
pred_obs_bwdo_fin$shelf <- factor(pred_obs_bwdo_fin$shelf, levels = c("West", "East"), labels = c("West", "East"))
summary(pred_obs_bwdo_fin$month)
pred_obs_bwdo_fin$monthf <- month.abb[pred_obs_bwdo_fin$month]
pred_obs_bwdo_fin$monthf <- parse_factor((pred_obs_bwdo_fin$monthf), levels = month.abb)
glimpse(pred_obs_bwdo_fin)
## Obs vs pred
(r2e_oh <- round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "East"], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "East"]), 2))
(r2w_oh <- round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "West"], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "West"]), 2))
(r2_oh <- round(E(pred_obs_bwdo_fin$Hindcasted, pred_obs_bwdo_fin$Observed), 2))
### Plot
ann_text_hosm <- data.frame(Observed= rep(2, 8), #1.5
                            Hindcasted = rep(5.2, 8), 
                            shelf = c(rep("West", 4), rep("East" , 4)),
                            monthf = c(rep(unique(pred_obs_bwdo_fin$monthf), 2)),
                            label = c(paste("italic(R) ^ 2 == ", round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "West" & pred_obs_bwdo_fin$month == 6], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "West" & pred_obs_bwdo_fin$month == 6]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "West" & pred_obs_bwdo_fin$month == 7], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "West" & pred_obs_bwdo_fin$month == 7]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "West" & pred_obs_bwdo_fin$month == 8], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "West" & pred_obs_bwdo_fin$month == 8]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "West" & pred_obs_bwdo_fin$month == 9], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "West" & pred_obs_bwdo_fin$month == 9]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "East" & pred_obs_bwdo_fin$month == 6], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "East" & pred_obs_bwdo_fin$month == 6]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "East" & pred_obs_bwdo_fin$month == 7], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "East" & pred_obs_bwdo_fin$month == 7]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "East" & pred_obs_bwdo_fin$month == 8], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "East" & pred_obs_bwdo_fin$month == 8]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(pred_obs_bwdo_fin$Hindcasted[pred_obs_bwdo_fin$shelf == "East" & pred_obs_bwdo_fin$month == 9], pred_obs_bwdo_fin$Observed[pred_obs_bwdo_fin$shelf == "East" & pred_obs_bwdo_fin$month == 9]), 2), sep = "")
                            ))
ann_text_hosm$monthf <- month.abb[ann_text_hosm$monthf]
ann_text_hosm$monthf <- parse_factor(ann_text_hosm$monthf, levels = month.abb)
ann_text_hosm$shelf <- factor(ann_text_hosm$shelf, levels = c("West", "East"))
glimpse(ann_text_hosm)
###
hind_obs_plot_two_shelves_monthly <- ggplot(pred_obs_bwdo_fin, aes(x = Hindcasted, Observed))+
  geom_point(shape = 21, alpha = 1)+
  facet_grid(shelf ~ monthf)+
  geom_abline(slope = 1, color = "grey50", size = 0.3)+
  theme_bw()+
  theme(legend.position = "top")+
  scale_color_manual(values = Palette_months())+
  coord_fixed(ratio = 1, xlim = c(1, 6), ylim = c(1, 6), expand = TRUE)+
  labs(color = "", y = expression(Observed~BWDO~(mg/L)), x = expression(Hindcasted~BWDO~(mg/L)))+
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
hind_obs_plot_two_shelves_monthly
##### only june bias correction
bias <- pred_obs_bwdo_fin %>% mutate(bias = (Observed - Hindcasted) / Hindcasted)
bias$monthday <- format(as.Date(bias$Date), "%m-%d")
bias$daynumber <- rep(rep(seq(1,122), 32), 2)
p_o_bwdo_june_west <- bias %>% filter(shelf == "West", month == 6)
glimpse(p_o_bwdo_june_west)
bias_plot <- p_o_bwdo_june_west %>% #filter(shelf == "West") %>% 
  ggplot(aes(daynumber, bias))+
  geom_point(shape = 1)+
  #facet_wrap(~ shelf)+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", se = F, color = "red", formula= y ~ 0 + x)+
  scale_x_continuous(breaks = c(1, 10, 20, 30))+
  #labs(x = "Day number", y = expression(frac(Observed-Predicted,Predicted)))+
  #theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5, hjust = 1))
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))+
  labs(x = "Day number", y = "Adjustment factor")
bias_plot
ggsave(bias_plot, file = "./graphics for manuscript/bias_regression.png", dpi = 500,
       width = 10, height = 7, units = "cm")
mod_bias <- (lm(bias ~ daynumber, data = p_o_bwdo_june_west))
summary(mod_bias)
saveRDS(mod_bias, "./data/bias_regression.rds")
bias_pred <- data.frame(bias_corr = predict(mod_bias, data.frame(daynumber = seq(1:30))), 
                        monthday = bias$monthday[1:30]) %>% 
  mutate(bias_corr = bias_corr + 1)
predict(mod_bias, newdata = data.frame(daynumber = c(2,20)), se.fit = T)
simulate(object = mod_bias, nsim = 1, newdata = data.frame(daynumber = c(30)))
(se <- sqrt(diag(vcov(mod_bias))))
predict(mod_bias, data.frame(daynumber = 2), interval = "prediction")
predict(mod_bias, data.frame(daynumber = 2), interval = "confidence")
### test 
N <- 100000
sims <- vector(mode = "logical", length = N)
sims <- (rep(coef(mod_bias)[1], N) + rnorm(n = N, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[1])) +
        (rep(coef(mod_bias)[2], N) + rnorm(n = N, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[2])) * 2
quantile(sims, probs = c(0.05, 0.95))
predict(mod_bias, data.frame(daynumber = 2), interval = "confidence", level = 0.90)
###
## Save bias adjustment factor 
saveRDS(bias_pred, "./data/bias_correction.rds")
############################################################################################
### Bias adjustment new
##### only june bias correction
bias <- pred_obs_bwdo_fin %>% mutate(bias = (Observed - Hindcasted) / Hindcasted)
bias$monthday <- format(as.Date(bias$Date), "%m-%d")
bias$daynumber <- rep(rep(seq(1,122), 32), 2)
backward_seq <- seq(from = 29, to = 0)
p_o_bwdo_june_west <- bias %>% filter(shelf == "West", month == 6) %>%
  mutate(daynumber_rev = rep(backward_seq, 32))
glimpse(p_o_bwdo_june_west)

bias_plot <- p_o_bwdo_june_west %>% #filter(shelf == "West") %>% 
  ggplot(aes(rev(daynumber_rev), bias))+
  geom_point(shape = 1)+
  #facet_wrap(~ shelf)+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", se = F, color = "red", formula= y ~ 0 + rev(x), fullrange = T)+
  scale_x_continuous(breaks = backward_seq[c(25, 15, 5)], labels = rev(backward_seq[c(25, 15, 5)]), 
                     expand = c(0.0, 0.0))+
  #labs(x = "Day number", y = expression(frac(Observed-Predicted,Predicted)))+
  #theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5, hjust = 1))
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))+
  labs(x = "Day number", y = "Adjustment factor")
bias_plot
############################################################################################
#ggsave(bias_plot, file = "./graphics for manuscript/bias_regression_no_intercept.png", dpi = 500,
#       width = 10, height = 7, units = "cm")
############################################################################################
mod_bias <- (lm(bias ~ 0 + daynumber_rev, data = p_o_bwdo_june_west))
summary(mod_bias)
saveRDS(mod_bias, "./data/bias_regression_new.rds")
(bias_pred <- data.frame(bias_corr = predict(mod_bias, data.frame(daynumber_rev = backward_seq)), 
                        monthday = bias$monthday[1:30]) %>% 
  mutate(bias_corr = bias_corr + 1))
## R-squared of this regression
E(pred = mod_bias$fitted.values, obs = mod_bias$model$bias)
predict(mod_bias, newdata = data.frame(daynumber_rev = c(2,20)), se.fit = T)
simulate(object = mod_bias, nsim = 1, newdata = data.frame(daynumber_rev = c(30)))
(se <- sqrt(diag(vcov(mod_bias))))
predict(mod_bias, data.frame(daynumber_rev = 2), interval = "prediction")
predict(mod_bias, data.frame(daynumber_rev = 2), interval = "confidence")
### test 
N <- 100000
sims <- vector(mode = "logical", length = N)
sims <- (rep(coef(mod_bias)[1], N) + rnorm(n = N, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[1])) * 2
quantile(sims, probs = c(0.05, 0.95))
predict(mod_bias, data.frame(daynumber_rev = 2), interval = "confidence", level = 0.90)
###
## Save bias adjustment factor 
saveRDS(bias_pred, "./data/bias_correction_new.rds")
############################################################################################
pred_obs_bwdo_fin$monthday <- format(as.Date(pred_obs_bwdo_fin$Date), "%m-%d")
out_finale <- left_join(pred_obs_bwdo_fin, bias_pred, by = "monthday")
glimpse(out_finale)
out_finale <- out_finale %>%
  mutate(bias_corr = if_else(is.na(bias_corr) == TRUE, 1, bias_corr),
         #Hindcasted = Hindcasted * bias_corr)
         Hindcasted = ifelse(shelf == "West", Hindcasted * bias_corr, Hindcasted))
### Plot updated
ann_text_hosm_upd <- data.frame(Observed= rep(1.5, 8),
                            Hindcasted = rep(5, 8),#5.2 
                            shelf = c(rep("West", 4), rep("East" , 4)),
                            monthf = c(rep(unique(out_finale$monthf), 2)),
                            label = c(paste("italic(R) ^ 2 == ", round(E(out_finale$Hindcasted[out_finale$shelf == "West" & out_finale$month == 6], out_finale$Observed[out_finale$shelf == "West" & out_finale$month == 6]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(out_finale$Hindcasted[out_finale$shelf == "West" & out_finale$month == 7], out_finale$Observed[out_finale$shelf == "West" & out_finale$month == 7]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(out_finale$Hindcasted[out_finale$shelf == "West" & out_finale$month == 8], out_finale$Observed[out_finale$shelf == "West" & out_finale$month == 8]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(out_finale$Hindcasted[out_finale$shelf == "West" & out_finale$month == 9], out_finale$Observed[out_finale$shelf == "West" & out_finale$month == 9]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(out_finale$Hindcasted[out_finale$shelf == "East" & out_finale$month == 6], out_finale$Observed[out_finale$shelf == "East" & out_finale$month == 6]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(out_finale$Hindcasted[out_finale$shelf == "East" & out_finale$month == 7], out_finale$Observed[out_finale$shelf == "East" & out_finale$month == 7]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(out_finale$Hindcasted[out_finale$shelf == "East" & out_finale$month == 8], out_finale$Observed[out_finale$shelf == "East" & out_finale$month == 8]), 2), sep = ""),
                                      paste("italic(R) ^ 2 == ", round(E(out_finale$Hindcasted[out_finale$shelf == "East" & out_finale$month == 9], out_finale$Observed[out_finale$shelf == "East" & out_finale$month == 9]), 2), sep = "")
                            ))
ann_text_hosm_upd$monthf <- month.abb[ann_text_hosm_upd$monthf]
ann_text_hosm_upd$monthf <- parse_factor(ann_text_hosm_upd$monthf, levels = month.abb)
ann_text_hosm_upd$shelf <- factor(ann_text_hosm_upd$shelf, levels = c("West", "East"))
glimpse(ann_text_hosm_upd)
###
hind_obs_plot_two_shelves_monthly_upd <- ggplot(out_finale, aes(x = Hindcasted, Observed))+
  geom_point(shape = 21, alpha = 1)+
  facet_grid(shelf ~ monthf)+
  geom_abline(slope = 1, color = "grey50", size = 0.3)+
  theme_bw()+
  theme(legend.position = "top")+
  scale_color_manual(values = Palette_months())+
  coord_fixed(ratio = 1, xlim = c(1, 6), ylim = c(1, 6), expand = TRUE)+
  labs(color = "", y = expression(Observed~BWDO~(mg/L)), x = expression(Hindcasted~BWDO~(mg/L)))+
  geom_text(data = ann_text_hosm_upd, label = ann_text_hosm_upd$label, parse = T, color = "black", size = 3)+
  theme(strip.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -10, -10))
hind_obs_plot_two_shelves_monthly_upd
########### COmbined

hind_obs_plot_two_shelves_monthly_both <- ggplot(pred_obs_bwdo_fin, aes(x = Hindcasted, Observed))+
  geom_point(shape = 21, alpha = 1)+
  facet_grid(shelf ~ monthf)+
  geom_abline(slope = 1, color = "grey50", size = 0.3)+
  theme_bw()+
  theme(legend.position = "top")+
  scale_color_manual(values = Palette_months())+
  coord_fixed(ratio = 1, xlim = c(1, 6), ylim = c(1, 6), expand = TRUE)+
  labs(color = "", y = expression(Observed~BWDO~(mg/L)), x = expression(Hindcasted~BWDO~(mg/L)))+
  geom_point(data = out_finale %>% filter(month == 6, shelf == "West"), aes(x = Hindcasted, Observed), shape = 21, alpha = 1, color = "red")+
  geom_text(data = ann_text_hosm_upd %>% filter(monthf == "Jun", shelf == "West"), label = ann_text_hosm_upd$label[1], parse = T, color = "red", size = 3)+
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
hind_obs_plot_two_shelves_monthly_both
ggsave(hind_obs_plot_two_shelves_monthly_both, file = "./graphics for manuscript/bias_adj_bwdo.png", dpi = 500,
       width = 18, height = 12, 
       units = "cm")
###########
glimpse(pred_bwdo)
load(file = "./data/DOvsArea_se18.RData")
pred_bwdo$monthday <- format(as.Date(pred_bwdo$Date), "%m-%d")
pred_bwdo <- left_join(pred_bwdo, bias_pred, by = "monthday") %>% mutate(month = month(Date))

pred_ha_shelves <- pred_bwdo %>% 
  mutate(ha_pred_e = lmHO.E$coef[1] + lmHO.E$coef[2] * bwdoE + lmHO.E$coef[3] * bwdoE ^ 2,
         ha_pred_w = lmHO.W$coef[1] + lmHO.W$coef[2] * bwdoW + lmHO.W$coef[3] * bwdoW ^ 2,
         ha_pred_w_ba = if_else(month == 6, 
                                lmHO.W$coef[1] + lmHO.W$coef[2] * bwdoW * bias_corr + lmHO.W$coef[3] * (bwdoW * bias_corr) ^ 2,
                                lmHO.W$coef[1] + lmHO.W$coef[2] * bwdoW + lmHO.W$coef[3] * bwdoW ^ 2)) 
ggplot(pred_ha_shelves, aes(ha_pred_w, ha_pred_w_ba))+geom_point()+facet_wrap(~ month) 
obs_ha_shelves <- Out.dat %>% dplyr::select(Date = 1, ha_obs_e = 7, ha_obs_w = 12)
out_plot_fhoE <- left_join(pred_ha_shelves %>% dplyr::select(Date, Hindcasted = ha_pred_e), 
                           obs_ha_shelves %>% dplyr::select(Date, Observed = ha_obs_e), by = "Date") %>% 
  mutate(shelf = "East")
out_plot_fhoW <- left_join(pred_ha_shelves %>% dplyr::select(Date, Hindcasted = ha_pred_w), 
                           obs_ha_shelves %>% dplyr::select(Date, Observed = ha_obs_w), by = "Date") %>%
  mutate(shelf = "West")
out_plot_fhoW_ba <- left_join(pred_ha_shelves %>% dplyr::select(Date, Hindcasted = ha_pred_w_ba), 
                           obs_ha_shelves %>% dplyr::select(Date, Observed = ha_obs_w), by = "Date") %>%
  mutate(shelf = "West") 
out_plot_fho_shelves <- bind_rows(out_plot_fhoE, out_plot_fhoW)
out_plot_fho_shelves_ba <- bind_rows(out_plot_fhoE, out_plot_fhoW_ba)
###################################################################################################
## Combined hindcast versus observed, by month
ho_shelves <- out_plot_fho_shelves %>% mutate(month = month(Date)) 
ho_shelves$monthday <- format(as.Date(ho_shelves$Date), "%m-%d")
ho_shelves$monthf <- month.abb[ho_shelves$month]
ho_shelves$monthf <- parse_factor((ho_shelves$monthf), levels = month.abb)
ho_shelves$shelf <- factor(ho_shelves$shelf, 
                           levels = c("West", "East"), 
                           ordered = TRUE)
ann_text_ho <- data.frame(Observed= rep(5000, 8),
                          Hindcasted = rep(22000, 8), 
                          shelf = c(rep("West", 4), rep("East" , 4)),
                          monthf = c(rep(unique(ho_shelves$monthf), 2)),
                          label = c(paste("italic(R) ^ 2 == ", round(E(ho_shelves$Hindcasted[ho_shelves$shelf == "West" & ho_shelves$month == 6], ho_shelves$Observed[ho_shelves$shelf == "West" & ho_shelves$month == 6]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves$Hindcasted[ho_shelves$shelf == "West" & ho_shelves$month == 7], ho_shelves$Observed[ho_shelves$shelf == "West" & ho_shelves$month == 7]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves$Hindcasted[ho_shelves$shelf == "West" & ho_shelves$month == 8], ho_shelves$Observed[ho_shelves$shelf == "West" & ho_shelves$month == 8]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves$Hindcasted[ho_shelves$shelf == "West" & ho_shelves$month == 9], ho_shelves$Observed[ho_shelves$shelf == "West" & ho_shelves$month == 9]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves$Hindcasted[ho_shelves$shelf == "East" & ho_shelves$month == 6], ho_shelves$Observed[ho_shelves$shelf == "East" & ho_shelves$month == 6]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves$Hindcasted[ho_shelves$shelf == "East" & ho_shelves$month == 7], ho_shelves$Observed[ho_shelves$shelf == "East" & ho_shelves$month == 7]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves$Hindcasted[ho_shelves$shelf == "East" & ho_shelves$month == 8], ho_shelves$Observed[ho_shelves$shelf == "East" & ho_shelves$month == 8]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves$Hindcasted[ho_shelves$shelf == "East" & ho_shelves$month == 9], ho_shelves$Observed[ho_shelves$shelf == "East" & ho_shelves$month == 9]), 2), sep = "")
                          ))
ann_text_ho$monthf <- month.abb[ann_text_ho$monthf]
ann_text_ho$monthf <- parse_factor(ann_text_ho$monthf, levels = month.abb)
ann_text_ho$shelf <- factor(ann_text_ho$shelf, levels = c("West", "East"), ordered = TRUE)
hind_obs_plot_two_shelves_ha_monthly <- ggplot(ho_shelves, aes(x = Hindcasted, Observed))+
  geom_point(shape = 21, alpha = 1)+
  facet_grid(shelf ~ monthf)+
  geom_abline(slope = 1, color = "grey50", size = 0.3)+
  theme_bw()+
  theme(legend.position = "top")+
  #scale_color_manual(values = Palette_months())+
  coord_fixed(ratio = 1, xlim = c(0, 27000), ylim = c(0, 27000), expand = TRUE)+
  labs(color = "", y = expression(Observed~HA~(km^{2})), x = expression(Hindcasted~HA~(km^{2})))+
  geom_text(data = ann_text_ho, label = ann_text_ho$label, parse = T, color = "black", size = 3)+
  theme(strip.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -10, -10))
hind_obs_plot_two_shelves_ha_monthly
### 
### Combined hindcast versus observed, by month BIAS ADJUSTED
ho_shelves_upd <- out_plot_fho_shelves_ba %>% mutate(month = month(Date)) 
ho_shelves_upd$monthday <- format(as.Date(ho_shelves_upd$Date), "%m-%d")
ho_shelves_upd$monthf <- month.abb[ho_shelves_upd$month]
ho_shelves_upd$monthf <- parse_factor((ho_shelves_upd$monthf), levels = month.abb)
ho_shelves_upd$shelf <- factor(ho_shelves_upd$shelf, 
                           levels = c("West", "East"), 
                           ordered = TRUE)
ann_text_ho_upd <- data.frame(Observed= rep(2000, 8),
                          Hindcasted = rep(21100, 8), 
                          shelf = c(rep("West", 4), rep("East" , 4)),
                          monthf = c(rep(unique(ho_shelves_upd$monthf), 2)),
                          label = c(paste("italic(R) ^ 2 == ", round(E(ho_shelves_upd$Hindcasted[ho_shelves_upd$shelf == "West" & ho_shelves_upd$month == 6], ho_shelves_upd$Observed[ho_shelves_upd$shelf == "West" & ho_shelves_upd$month == 6]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves_upd$Hindcasted[ho_shelves_upd$shelf == "West" & ho_shelves_upd$month == 7], ho_shelves_upd$Observed[ho_shelves_upd$shelf == "West" & ho_shelves_upd$month == 7]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves_upd$Hindcasted[ho_shelves_upd$shelf == "West" & ho_shelves_upd$month == 8], ho_shelves_upd$Observed[ho_shelves_upd$shelf == "West" & ho_shelves_upd$month == 8]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves_upd$Hindcasted[ho_shelves_upd$shelf == "West" & ho_shelves_upd$month == 9], ho_shelves_upd$Observed[ho_shelves_upd$shelf == "West" & ho_shelves_upd$month == 9]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves_upd$Hindcasted[ho_shelves_upd$shelf == "East" & ho_shelves_upd$month == 6], ho_shelves_upd$Observed[ho_shelves_upd$shelf == "East" & ho_shelves_upd$month == 6]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves_upd$Hindcasted[ho_shelves_upd$shelf == "East" & ho_shelves_upd$month == 7], ho_shelves_upd$Observed[ho_shelves_upd$shelf == "East" & ho_shelves_upd$month == 7]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves_upd$Hindcasted[ho_shelves_upd$shelf == "East" & ho_shelves_upd$month == 8], ho_shelves_upd$Observed[ho_shelves_upd$shelf == "East" & ho_shelves_upd$month == 8]), 2), sep = ""),
                                    paste("italic(R) ^ 2 == ", round(E(ho_shelves_upd$Hindcasted[ho_shelves_upd$shelf == "East" & ho_shelves_upd$month == 9], ho_shelves_upd$Observed[ho_shelves_upd$shelf == "East" & ho_shelves_upd$month == 9]), 2), sep = "")
                          ))
ann_text_ho_upd$monthf <- month.abb[ann_text_ho_upd$monthf]
ann_text_ho_upd$monthf <- parse_factor(ann_text_ho_upd$monthf, levels = month.abb)
ann_text_ho_upd$shelf <- factor(ann_text_ho_upd$shelf, levels = c("West", "East"), ordered = TRUE)
hind_obs_plot_two_shelves_ha_monthly_ba <- ggplot(ho_shelves_upd, aes(x = Hindcasted, Observed))+
  geom_point(shape = 21, alpha = 1)+
  facet_grid(shelf ~ monthf)+
  geom_abline(slope = 1, color = "grey50", size = 0.3)+
  theme_bw()+
  theme(legend.position = "top")+
  #scale_color_manual(values = Palette_months())+
  coord_fixed(ratio = 1, xlim = c(0, 27000), ylim = c(0, 27000), expand = TRUE)+
  labs(color = "", y = expression(Observed~HA~(km^{2})), x = expression(Hindcasted~HA~(km^{2})))+
  geom_text(data = ann_text_ho_upd, label = ann_text_ho_upd$label, parse = T, color = "black", size = 3)+
  theme(strip.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -10, -10))
hind_obs_plot_two_shelves_ha_monthly_ba
######################################################
hind_obs_plot_two_shelves_ha_monthly_both <- ggplot(ho_shelves, aes(x = Hindcasted/10000, Observed/10000))+
  geom_point(shape = 21, alpha = 1)+
  geom_point(data = ho_shelves_upd %>% filter(month == 6, shelf == "West"), aes(x = Hindcasted/10000, Observed/10000), shape = 21, alpha = 1, color = "red")+
  facet_grid(shelf ~ monthf)+
  geom_abline(slope = 1, color = "grey50", size = 0.3)+
  theme_bw()+
  theme(legend.position = "top")+
  #scale_color_manual(values = Palette_months())+
  coord_fixed(ratio = 1, xlim = c(0, 2.7), ylim = c(0, 2.7), expand = TRUE)+
  labs(color = "", y = expression(Observed~HA~(10^{4}*km^{2})), x = expression(Hindcasted~HA~(10^{4}*km^{2})))+
  geom_text(data = ann_text_ho, label = ann_text_ho$label, parse = T, color = "black", size = 3)+
  geom_text(data = ann_text_ho_upd %>% filter(monthf == "Jun", shelf == "West"), label = ann_text_ho_upd$label[1], parse = T, color = "red", size = 3)+
  theme(strip.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -10, -10))
hind_obs_plot_two_shelves_ha_monthly_both
ggsave(hind_obs_plot_two_shelves_ha_monthly_both, file = "./graphics for manuscript/bias_adj_ha.png", dpi = 500,
       width = 18, height = 12, 
       units = "cm")
###################################################################