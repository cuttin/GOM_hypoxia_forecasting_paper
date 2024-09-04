if (!require("pacman")) install.packages("pacman"); 
pacman::p_load(magrittr, tidyverse, lubridate, MHadaptive, IDPmisc, reshape2, tictoc, ggthemes, 
               faraway, foreach, doParallel, rstudioapi) 
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
wT0 <- read_csv(file="shelf_salinity_temp_apr18.csv", skip=0)
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
### Load model parameters:
#load("./data/GmCalbT_ap19.RData")
#mcmc <- get("MHa")
#str(mcmc)
#par_all <- (mcmc$trace[seq((nrow(mcmc$trace)-2002), nrow(mcmc$trace), 2),]) # best parameter set
#saveRDS(par_all, "./data/pars.rds")
par_all <- readRDS(file = "./data/pars.rds")
par_mean <- colMeans(par_all)
# fixed parameters, hardcoded in the model 
R_CN = 5.68; R_OC = 3.5; # ratio of carbon/nitr, oxyg/carbon
A.W = 48450 / 1000  #shelf area 1000km2, G(m2)
A.E = 14000 / 1000  #shelf area 1000km2, G(m2)
CS_O = 3 #reference DO for sediment oxygen demand
FML = 0.2 #fraction of MS River lost under all wind conditions
DOadj = 1 #adjustment to model lower layer DO, not BWDO
rQg = 3.2 # Ocean dilution factor
meanMsQ = mean(colMeans(spgMsQ))
###############################################################################################
############ Function to run the model
############ Vector one year run, this is way faster compared to code in DMO
one_year_run_vector <- function(data, params) {
  data <- data %>% add_column(bwdoE = NA, bwdoW = NA, bwdoW_e = NA, bwdoE_e = NA)
  
  # East/Mississ bottom segm
  pngFe    <- 0.5 + params["bE"] * data$spgWwi #eq 5, spring / WCOD related 
  pngFe    <- pmin(pmax(1e-2, pngFe), 1) # Fe constrained
  smrFe    <- 0.5 + params["bE"] * data$smrWwi #eq 5, summer / reaeration related
  smrFe    <- pmin(pmax(1e-2, smrFe), 1)
  QsE      <- (data$smrAtQ * 86400 / 10 ^ 9) * smrFe # for eq 6, summer freshwater input [Gm3/d]
  tauE     <- (data$smr14wi2W / 360 + data$smr14wi2E / 110) / (1 / 360 + 1 / 110) / 100
  kaE      <- params["b0Ka"] + params["b1Ka"] * tauE / QsE * A.E / 100 #eq. 6, summer reaerat [m/d]
  NinE     <- data$spgMsN / 30.5 / 1000 * (1 - FML) * (1 - pngFe) + data$spgAtN / 30.5 / 1000 * pngFe # load entering in spring (QrCrN+QuCuN) [Gg/d]
  QinE     <- (rQg * meanMsQ * 86400 / 10 ^ 9 + 
                 data$spgMsQ * 86400 / 10 ^ 9 * (1 - FML)) * (1 - pngFe) + 
    data$spgAtQ * 86400 / 10 ^ 9 * pngFe # [Gm3/d]
  CfxE  <- NinE * R_CN * params["vs"] / (A.E * params["vs"] + QinE) # gC/(m2*d)
  data$bwdoE <- (kaE * data$Os - CfxE * R_OC * params["R_PR"]) / (kaE + (data$LoNoMa / avgNwint) ^ (.5) *
                                                                    params["SD"] * 1.07 ^ (data$bWt - ref_bwT) / CS_O) - DOadj
  
  # West/Atchafa bottom segm
  QsW      <- (data$smrAtQ * 86400 / 10 ^ 9) * (1 - smrFe) # for eq 6, [Gm3/d]
  tauW     <- (data$smr14wi2W / 130 + data$smr14wi2E / 350)/ (1 / 130 + 1 / 350) / 100
  kaW      <- params["b0Ka"] + params["b1Ka"] * tauW / QsW * A.W / 100 # eq. 6, summer reaerat [m/d]
  NoutE    <- NinE / (A.E * params["vs"] + QinE) * QinE # reduced by setting [Gg/d]
  NinW     <- (NoutE + data$spgAtN / 30.5 / 1000) * (1 - pngFe) # available loads moved westward
  QinW     <- (QinE + data$spgAtQ * 86400 / 10^9) * (1 - pngFe) + 
    rQg * meanMsQ * 86400 / 10 ^ 9 * pngFe # [Gm3/d]
  CfxW     <- NinW * R_CN * params["vs"] / (A.W * params["vs"] + QinW)
  data$bwdoW <- (kaW * data$Os - CfxW * R_OC * params["R_PR"]) / (kaW + (data$LoNoMa / avgNwint) ^ (.5) * 
                                                                    params["SD"] * 1.07 ^ (data$bWt - ref_bwT) / CS_O) - DOadj #
  
  data$bwdoW_e <- data$bwdoW + rnorm(nrow(data), mean = 0, sd = params["sy.W"])
  data$bwdoE_e <- data$bwdoE + rnorm(nrow(data), mean = 0, sd = params["sy.E"])
  data <- data %>% dplyr::select(Date, bwdoW, bwdoW_e, bwdoE, bwdoE_e)
  return(data)
}
## check
tic()
dd1 <- data.frame()
for (i in seq_along(years_of_calibration)) {
  tem <- one_year_run_vector(data = InpDaif %>% filter(year == years_of_calibration[i]), params = par_mean)
  dd1 <- bind_rows(dd1, tem)
}
toc()

apply(par_all[1:50, ], 1, function(x) one_year_run_vector(data = InpDaif %>% filter(year == years_of_calibration[1]), params = x))

#####################################################################################
RivDai <- read_csv(file = "./data/daiAtMs_feb18.csv")
RivDai$time <- seq(as.Date("1980/01/1"), as.Date("2016/09/30"), "days")
RivDai <- RivDai %>% rename(Date = time) 
RivDai[ , c(3, 5)] <- RivDai[ , c(3, 5)] * 30.5 * 86400 # new conversion to /mo
respnames <- c("AtQ", "AtN", "MsQ", "MsN")
colnames(RivDai)[c(2, 3, 4, 5)] <- respnames
glimpse(RivDai)
## Read daily winds
wind_dates <- seq.Date(from = as.Date("1985/01/01"), to = as.Date("2016/12/31"), by = "1 day")
winds_ee <- read_csv(file = "./data/dayly_winds_e.csv", col_names = F)
winds_ww <- read_csv(file = "./data/dayly_winds_w.csv", col_names = F)
winds_e <- bind_cols(Date = wind_dates, winds_ee) %>% 
  dplyr::select(Date, year = X2, month = X3, ece = X4, nce = X5, wse = X6, missinge = X7)
winds_w <- bind_cols(Date = wind_dates, winds_ww) %>% 
  dplyr::select(Date, ecw = X4, ncw = X5, wsw = X6, missingw = X7)
winds <- left_join(winds_e, winds_w, by = "Date") %>% 
  filter(month %in% c(3:9))
winds$monthday <- factor(format(winds$Date, format = "%m-%d"))
glimpse(winds)
# weights for 14-day windspeeds 
weights_winds <- (1:14/(sum(1:14)))
## Read interpolated monthly winds
intp_e <- read_csv('./data/intp_e.csv')
intp_w <- read_csv('./data/intp_w.csv')
glimpse(intp_w)
intp_winds <- left_join(intp_e, intp_w, by = "Date")
## Prepare for input
qn_inp <- RivDai %>% mutate(year = year(Date),
                            month = month(Date)) %>% 
  filter(year >= 1985) %>% # filter year > 1985
  filter(month %in% c(3:9))
## 
min_daily <- qn_inp %>% 
  filter(month %in% c(6:9)) %>%
  group_by(month) %>%
  summarise_all(min) %>%
  dplyr::select(-Date, -year) %>%
  melt(., id = "month") %>%
  reshape2::dcast(., . ~ variable + month , value.var = "value")
## Monthly values
RivDaim <- RivDai %>% 
  mutate(year = year(Date),
         month = month(Date)) %>%
  filter(month %in% c(4:9)) %>%
  filter(year >= 1985) %>% # filter year > 1985
  dplyr::select(-Date) %>% 
  group_by(year, month) %>% 
  summarise_all(mean) %>%
  melt(., id = c("year", "month"))
##
dailyseq <- seq(as.Date("1967/10/1"), as.Date("2016/9/1"), "day")
monthlyseq <- data.frame(Date = dailyseq) %>% 
  mutate(year = year(Date), month = month(Date)) %>%
  group_by(year, month) %>%
  summarise(daysinmonth = as.numeric(n())) # also counting number of days in a month
gommonthly <- readxl::read_xlsx("./data/Gulf-Monthly-2016.xlsx", skip = 40, col_names = F)
gomm <- bind_cols(monthlyseq, gommonthly %>% dplyr::select(AtQm = 6, AtNOx = 57, AtTKN = 61, AtNH3 = 65, 
                                                    MsQm = 3, MsNOx = 32, MsTKN = 37, MsNH3 =  41))
gom <- gomm %>% 
  mutate(AtNbio = (AtNOx + AtNH3 + (AtTKN - AtNH3) * .12), # tonns / second
         MsNbio = (MsNOx + MsNH3 + (MsTKN - MsNH3) * .12),
         AtCm = AtNbio / AtQm,
         MsCm = MsNbio / MsQm) %>% # this is in mega gramms (gamms * 10^-6) because loading was in metric tonns
  dplyr::select(year, month, AtQ = AtQm, AtN = AtNbio, MsQ = MsQm, MsN = MsNbio) %>%
  filter(year > 1979)
ggplot(gom, aes(AtN/AtQ, MsN/MsQ))+geom_point()
head(gom)
glimpse(gom)
RivDaim <- gom
RivDaim <- melt(RivDaim, id = c("year", "month"))
glimpse(RivDaim)
## Without 1:4
RivDaimm <- gom %>%
  filter(month %in% c(5:9)) %>%
  melt(., id = c("year", "month"))
glimpse(RivDaimm)
jfma <- RivDaim %>% 
  filter(month %in% c(1, 2, 3, 4)) %>% 
  group_by(year, variable) %>%
  summarise_all(mean) %>%
  dplyr::select(colnames(RivDaim))
jfma$month <- if_else(jfma$month == 2.5, 1234, jfma$month)
reg_input <- bind_rows(RivDaimm, jfma) %>% 
  arrange(year, variable, month) %>%
  reshape2::dcast(., year ~ variable + month , value.var = "value")
## Read precipitation and temperature
precip <- read_csv("./data/prism_prec_monthly.csv")
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
predictors <- left_join(win, may, by = "year")
glimpse(predictors)
reg_input <- left_join(reg_input, predictors, by = "year")
glimpse(reg_input)
## load_models
load("./data/AMQL_models_sqrt_prism_05-11.RData")
pr <- data.frame(preds = matrix(NA, nrow = 16, ncol = 1))
for (i in seq_along(mod_list)) {
  pr[i, ] <- paste0(unlist(mod_list[[i]][[5]]), collapse = ", ")
}
tibble(response = sapply(mod_list, function(x) x[[1]]),
       predictors = pr$preds,
       adj.rsquared = round(sapply(mod_list, function(x) x[[4]]), 2)
)
#sel_list <- mod_list[c(1, 2, 5, 6, 9, 10, 13, 14)] # select models for june and july
sel_list <- mod_list
sel_list[[4]]
#lm_resp <- c("AtQ_6", "AtQ_7", "AtN_6", "AtN_7", "MsQ_6", "MsQ_7", "MsN_6", "MsN_7" )
lm_resp <- sapply(mod_list, function(x) x[[1]])
###
# set years for predictions 
years_of_interest <- 1995
years_of_interest <- years_of_calibration
#years_of_interest <- 1998:1999
# how Many simulations 
N = 1000
## Create vector with list 
list_years <- list()
list_pred_vals <- list()
list_fin <- list()
list_fin_hm <- list() # save all BWDO sims only meteorology
#hydrometeo_all <- list()
list_sims_all <- list()
negsims <- matrix(NA, nrow = 16, ncol = length(years_of_interest))
#################################################################################################
#################################################################################################
#################################################################################################
######### 10 relevant years
for (j in seq_along(years_of_interest)) {
  #for (j in seq_along(years_of_interest[years_of_interest == 2005])) {
  yy <- years_of_interest[j]
  yyno <- years_of_calibration[years_of_calibration != yy]
  inp <- reg_input %>% filter(year == yy) 
    # dataframe with year to model and only spring data 
  inpout <- reg_input %>% filter(year != yy) # dataframe with other years
  #  
  list_years_all <- list()
  list_sims <- list()
  predvals <- data.frame()
    for (k in 1:length(sel_list)) {
    
    cmod <- sel_list[[k]]
    varname <- cmod$response
    #mind <- min_daily %>% dplyr::pull(varname) # get minimum daily observed value for constraining predictions
    #datain <- inp[, colnames(inp) %in% names(cmod$model$coefficients)]
    #res <- data.frame(predict(cmod$model, inp, interval = "prediction", level = 0.95)) %>%
    #                  transmute(var = varname, !!!.) %>%
    #  mutate(fit = fit ^ 2, lwr = lwr, upr = upr)# predicted variable column comes first
    res <- data.frame(predict(cmod$model, inp, interval = "prediction", level = 0.95))
    res <- res %>% dplyr::mutate(var = varname, .before = fit) %>%
      mutate(fit = fit ^ 2, lwr = lwr, upr = upr)# predicted variable column comes first
    cmod$model$fitted.values <- predict(cmod$model, inp) # "hack"
    # simulated_values <- vector(length = N, mode = "numeric")
    # for (i in seq_along(simulated_values)) {
    #   repeat {
    #     simulated_values[i] <- as.numeric(simulate(cmod$model, nsim = 1, interval = "prediction"))
    #   if (simulated_values[i] >= mind){
    #     break
    #         }
    #     }
    # }
    
    simulated_values <- melt(simulate(cmod$model, nsim = N, seed = 123)) %>% 
      mutate(value = value ^ 2) %>%
      dplyr::pull(value)
    negsims[k, j] <- length(which(simulated_values <= 0))
    #simulated_values <- if_else(simulated_values <= 0, sample(simulated_values[simulated_values > 0], 1), simulated_values)
    #observed_value <- inp %>% pull(varname) # extract observed value
    #res <- mutate(res, f = observed_value / fit) # calculate fraction observed vs 
    vari <- inpout %>% dplyr::select(year, varname)  # get the observed values of response variable
    # get sd and mean of the data (all years)
    meanvar <- mean(reg_input %>% dplyr::pull(varname), na.rm = T) 
    sdvar <- sd(reg_input %>% dplyr::pull(varname), na.rm = T)
    #yyy <- vari$year[data.table::between(vari[, 2], res[, 3], res[, 4])] # see if value is between using interval
    # get difference between scaled and centered observed values for other years and 
    # predicted value for forecasted year 
    vari$dif <- abs((vari[, 2] - meanvar) / sdvar - (res[, 2] - meanvar) / sdvar) 
    yyy <- (vari %>% arrange(dif) %>% dplyr::select(year, dif) %>% mutate(var = varname))[1:10, ] # extract those years
    predvals <- bind_rows(predvals, res)
    list_years_all[[k]] <- yyy 
    list_sims[[k]] <- simulated_values
    names(list_sims)[k] <- varname # give a name to each set of simulations
        } 
  #tempy <- Reduce(rbind, list_years_all)
  tempy <- Reduce(rbind, list_years_all[c(1, 2, 5, 6, 9, 10, 13, 14)])
  tempy <- dplyr::arrange(tempy, dif)
  tempy <- tempy[tempy$year %in% names(which(table(tempy$year) > 2)), ] # remove years which are not there at least 2 times
  tempyd <- (dplyr::distinct(tempy, year, .keep_all = TRUE))[1:10, ]
  #list_years[[j]] <- Reduce(intersect, list_years_all) # select matching years
  #list_years[[j]] <- sort(unique(unlist(list_years_all[c(1, 2, 5, 6, 9, 10, 13, 14)]))) # select all unique years for June-July models 
  list_years[[j]] <- tempyd %>% 
    dplyr::arrange(year) %>% 
    dplyr::pull(year)
  list_pred_vals[[j]] <- predvals
  list_sims_all[[j]] <- list_sims 
  names(list_years)[j] <- years_of_interest[j]
  names(list_pred_vals)[j] <- years_of_interest[j]
  names(list_sims_all)[j] <- years_of_interest[j]
  print(years_of_interest[j])
}  

list_years
list_pred_vals
list_years_all
list_sims_all[[1]]
#################################################################################################
#################################################################################################
#################################################################################################
##################### Run data generation, obtain hydrometeorological inputs
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
foreach::foreach(j = seq_along(years_of_calibration), 
                 .packages = c('dplyr', 'lubridate', 'tibble')) %dopar% {
                 #.export = c('logpri', 'loglik', 'logpost', 
                  #           'Metro_Hastings_0815', 'prior.pbdis', 'par')) %dopar% {

  yy <- years_of_interest[j] # year to predict
  yyno <- years_of_calibration[years_of_calibration != yy] # all years except year to predict
  inp <- reg_input %>% dplyr::filter(year == yy) # dataframe with year to model
  inpout <- reg_input %>% dplyr::filter(year != yy) # dataframe with other years
  dates_to_model <- seq(from = as.Date(paste(yy, "06", "01",sep = "-", format = "%y-%m-%d")),
                        to = as.Date(paste(yy, "09", "30",sep = "-", format = "%y-%m-%d")),
                        by = "1 day")
  qn_actual <- dplyr::filter(qn_inp, year == yy)
  qn_actual_spring <- dplyr::filter(qn_actual, month < 6, year == yy)
  qn_actual_spring$monthday <- factor(format(qn_actual_spring$Date, format = "%m-%d"))
  # get actual winds
  ##########################################################
  #!!MUST BE < 6
  winds_actual_spring <- winds %>% dplyr::filter(month < 6, year == yy)
  #  
  const_inp <- InpDaif %>% 
    dplyr::filter(year == yy) %>% 
    dplyr::select(Date, LoNoMa, Os, bWt)
  
  
  # list of outputs
  hydrometeo <- list()
  
  for (x in seq_along(1:N)) {
    
    # create empty df
    InpDaiGomMod <- data.frame(matrix(NA, nrow = length(dates_to_model), ncol = 5))
    colnames(InpDaiGomMod) <- c("spgAtQ", "spgAtN", "spgMsQ", "spgMsN", "smrAtQ")
    InpDaiGomMod <- InpDaiGomMod %>% 
      tibble::add_column(spgWwi = NA, smrWwi = NA, smr14wi2W = NA, smr14wi2E = NA)
    InpDaiGomMod <- dplyr::bind_cols(Date = dates_to_model, InpDaiGomMod)
    InpDaiGomMod <- dplyr::bind_cols(InpDaiGomMod, const_inp %>% 
                                       dplyr::select(-Date))
    
    picked_year <- sample(list_years[[j]], 1, replace = T) # sample one year from the history
    fs <- reg_input %>% dplyr::filter(year == picked_year) # mean monthly values for the year
    other_qn_inp <- RivDai %>% 
      dplyr::mutate(year = year(Date), month = month(Date)) %>%
      dplyr::filter(year == picked_year, month >= 6) # leave only daily inputs for one of the prev years
    head(other_qn_inp)
    other_qn_inp_upd <- other_qn_inp 
    other_qn_inp_upd$monthday <- factor(format(other_qn_inp_upd$Date, format = "%m-%d"))
    
    random_year <- sample(yyno, 1, replace = T) # get a random year
    other_winds <- winds %>% dplyr::filter(month >= 6, year == random_year)
    
    # start loop for flows and loadings
    for (i in 1:length(dates_to_model) ) {
      
      back90 <- seq(dates_to_model[i], by = "-1 days", length.out = 91) 
      
      per1_90_60 <- back90[62:91] # used only for spring (2-3 months before pred point/cruise)
      inp_90_60 <- qn_actual_spring %>% dplyr::filter(Date %in% per1_90_60)  
      future_dates_90_60 <- data.frame(Date = per1_90_60[per1_90_60 %in% dates_to_model])
      future_dates_90_60$monthday <- factor(format(future_dates_90_60$Date, format = "%m-%d"))
      from_prev_years_90_60 <- other_qn_inp_upd %>% 
        dplyr::filter(monthday %in% future_dates_90_60$monthday) %>%
        dplyr::select(-monthday)
      inp_90_60 <- dplyr::bind_rows(inp_90_60, from_prev_years_90_60)   
      
      per2_60_30 <- back90[32:61] # used for both spring and summer
      
      inp_60_30 <- qn_actual_spring %>% dplyr::filter(Date %in% per2_60_30)  
      future_dates_60_30 <- data.frame(Date = per2_60_30[per2_60_30 %in% dates_to_model])
      future_dates_60_30$monthday <- factor(format(future_dates_60_30$Date, format = "%m-%d"))
      from_prev_years_60_30 <- other_qn_inp_upd %>% 
        dplyr::filter(monthday %in% future_dates_60_30$monthday) %>%
        dplyr::select(-monthday)
      inp_60_30 <- dplyr::bind_rows(inp_60_30, from_prev_years_60_30)
      
      per3_00_30 <- back90[1:31] # used only for summer (1-2 months before pred point/cruise)
      inp_00_30 <- qn_actual_spring %>% dplyr::filter(Date %in% per3_00_30)  
      future_dates_00_30 <- data.frame(Date = per3_00_30[per3_00_30 %in% dates_to_model])
      future_dates_00_30$monthday <- factor(format(future_dates_00_30$Date, format = "%m-%d"))
      from_prev_years_00_30 <- other_qn_inp_upd %>% 
        dplyr::filter(monthday %in% future_dates_00_30$monthday) %>%
        dplyr::select(-monthday)
      inp_00_30 <- dplyr::bind_rows(inp_00_30, from_prev_years_00_30)
      
      InpDaiGomMod$spgAtQ[i] <- (mean(inp_90_60$AtQ) * 2 + 
                                   mean(inp_60_30$AtQ)) / 3 # m3/s 
      
      InpDaiGomMod$spgAtN[i] <- (mean(inp_90_60$AtN) * 2 + 
                                   mean(inp_60_30$AtN)) / 3 # Mg/mo
      
      InpDaiGomMod$spgMsQ[i] <- (mean(inp_90_60$MsQ) * 2 + 
                                   mean(inp_60_30$MsQ)) / 3 # m3/s
      
      InpDaiGomMod$spgMsN[i] <- (mean(inp_90_60$MsN) * 2 + 
                                   mean(inp_60_30$MsN)) / 3 # Mg/mo
      
      InpDaiGomMod$smrAtQ[i] <- (mean(inp_00_30$AtQ) * 2 + 
                                   mean(inp_60_30$AtQ))/3 # m3/s
      
      ## Winds
      
      winds_90_60 <- winds_actual_spring %>% 
        dplyr::filter(Date %in% per1_90_60)
      winds_90_60 <- dplyr::left_join(winds_90_60, intp_winds, by = "Date") %>%
        dplyr::mutate(ecef = ece * missinge + (1 - missinge) * interp_ece,
               ecwf = ecw * missingw + (1 - missingw) * interp_ecw)
      future_dates_90_60 <- data.frame(Date = per1_90_60[per1_90_60 %in% dates_to_model])
      future_dates_90_60$monthday <- factor(format(future_dates_90_60$Date, format = "%m-%d"))
      winds_prev_years_90_60 <- other_winds %>% 
        dplyr::filter(monthday %in% future_dates_90_60$monthday)
      winds_prev_years_90_60 <- dplyr::left_join(winds_prev_years_90_60, intp_winds, by = "Date") %>%
        dplyr::mutate(ecef = ece * missinge + (1 - missinge) * interp_ece,
                      ecwf = ecw * missingw + (1 - missingw) * interp_ecw)
      winds_90_60 <- dplyr::bind_rows(winds_90_60, winds_prev_years_90_60)
      
      winds_60_30 <- winds_actual_spring %>% dplyr::filter(Date %in% per2_60_30)
      winds_60_30 <- dplyr::left_join(winds_60_30, intp_winds, by = "Date") %>%
        dplyr::mutate(ecef = ece * missinge + (1 - missinge) * interp_ece,
                      ecwf = ecw * missingw + (1 - missingw) * interp_ecw)
      future_dates_60_30 <- data.frame(Date = per2_60_30[per2_60_30 %in% dates_to_model])
      future_dates_60_30$monthday <- factor(format(future_dates_60_30$Date, format = "%m-%d"))
      winds_prev_years_60_30 <- other_winds %>% 
        dplyr::filter(monthday %in% future_dates_60_30$monthday)
      winds_prev_years_60_30 <- dplyr::left_join(winds_prev_years_60_30, intp_winds, by = "Date") %>%
        dplyr::mutate(ecef = ece * missinge + (1 - missinge) * interp_ece,
                      ecwf = ecw * missingw + (1 - missingw) * interp_ecw)
      winds_60_30 <- dplyr::bind_rows(winds_60_30, winds_prev_years_60_30)
      
      winds_00_30 <- winds_actual_spring %>% 
        dplyr::filter(Date %in% per3_00_30)
      winds_00_30 <- dplyr::left_join(winds_00_30, intp_winds, by = "Date") %>%
        dplyr::mutate(ecef = ece * missinge + (1 - missinge) * interp_ece,
                      ecwf = ecw * missingw + (1 - missingw) * interp_ecw)
      future_dates_00_30 <- data.frame(Date = per3_00_30[per3_00_30 %in% dates_to_model])
      future_dates_00_30$monthday <- factor(format(future_dates_00_30$Date, format = "%m-%d"))
      winds_prev_years_00_30 <- other_winds %>% dplyr::filter(monthday %in% future_dates_00_30$monthday)
      winds_prev_years_00_30 <- dplyr::left_join(winds_prev_years_00_30, intp_winds, by = "Date") %>%
        dplyr::mutate(ecef = ece * missinge + (1 - missinge) * interp_ece,
                      ecwf = ecw * missingw + (1 - missingw) * interp_ecw)
      winds_00_30 <- dplyr::bind_rows(winds_00_30, winds_prev_years_00_30)
      
      per4_00_14 <- back90[2:15]
      winds_00_14 <- winds_actual_spring %>% 
        dplyr::filter(Date %in% per4_00_14)
      winds_00_14 <- dplyr::left_join(winds_00_14, intp_winds, by = "Date") %>%
        dplyr::mutate(wsef = wse * missinge + (1 - missinge) * interp_wse,
                      wswf = wsw * missingw + (1 - missingw) * interp_wsw)
      future_dates_00_14 <- data.frame(Date = per4_00_14[per4_00_14 %in% dates_to_model])
      future_dates_00_14$monthday <- factor(format(future_dates_00_14$Date, format = "%m-%d"))
      winds_prev_years_00_14 <- other_winds %>% 
        dplyr::filter(monthday %in% future_dates_00_14$monthday)
      winds_prev_years_00_14 <- dplyr::left_join(winds_prev_years_00_14, intp_winds, by = "Date") %>%
        dplyr::mutate(wsef = wse * missinge + (1 - missinge) * interp_wse,
                      wswf = wsw * missingw + (1 - missingw) * interp_wsw)
      winds_00_14 <- dplyr::bind_rows(winds_00_14 %>% dplyr::select(-monthday), winds_prev_years_00_14)
      winds_00_14$weights_calc <- weights_winds
      winds_00_14 <- winds_00_14 %>% 
        dplyr::mutate(wei_wse = (wsef ^ 2) * weights_calc, wei_wsw = (wswf ^ 2) * weights_calc)
      
      InpDaiGomMod$spgWwi[i] <- (mean(append(winds_90_60$ecef, winds_90_60$ecwf)) * 2 + 
                            mean(append(winds_60_30$ecef, winds_60_30$ecwf))) / 3 # m3/s
      
      InpDaiGomMod$smrWwi[i] <- (mean(append(winds_00_30$ecef, winds_00_30$ecwf)) * 2 + 
                            mean(append(winds_60_30$ecef, winds_60_30$ecwf))) / 3 # m3/s
      
      InpDaiGomMod$smr14wi2W[i] <- sum(winds_00_14$wei_wsw)
      InpDaiGomMod$smr14wi2E[i] <- sum(winds_00_14$wei_wse)
      
    } # end of loop for daily flows loading and winds 
    
    hydrometeo[[x]] <- InpDaiGomMod
    
    } # closing loop for one simulation out of N
  saveRDS(hydrometeo, paste0("./data/hydrometeo_sims_01-21_paral/", yy, ".rds"))
  print(yy)
} # closing loop for the forecasting year

parallel::stopCluster(cl)

### test
yy <- 2009
test_hydro <- readRDS(file = paste0("./data/hydrometeo_sims_01-21_paral/", yy, ".rds"))
str(test_hydro[[1]])
id_col <- rep(1:1000, times = 1, each = 122)
test_hydro_dat <- do.call(bind_rows, test_hydro)
test_hydro_dat$id <- id_col
ggplot(test_hydro_dat %>% filter(id %in% c(1:100)), aes(Date, spgMsN, color = id))+
  geom_point(alpha = 0.2)
########################################################################################
########################################################################################
###
### Simulations
memory.limit()
#memory.limit(60000)
## Function to get mean and 2.5 and 97.5 quantiles
p95 <- c(0.025, 0.975)
p_names95 <- map_chr(p95, ~paste0(.x * 100, "%"))
p_funs95 <- map(p95, ~ partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names95)
## Load bias correction factor, no parameter uncertainty 
bias_pred <- readRDS(file = "./data/bias_correction.rds")
bias_vector <- vector(mode = 'numeric', length = 122L)
bias_vector[1:30] <- bias_pred[, 1]
bias_vector <- dplyr::if_else(bias_vector == 0, 1, bias_vector)
####
#### Get regression uncertainty
mod_bias <- readRDS("./data/bias_regression.rds")
summary(mod_bias)
bias_sims <- seq(1, 122)
(bias_sims <- ifelse(bias_sims > 30, NA, bias_sims))
N <- 1000
bias_1000_sims <- vector()
bias_1000_sims_list <- list()
for (i in 1:N) {
sims <- 1 + (coef(mod_bias)[1] + rnorm(n = 1, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[1])) +
        (coef(mod_bias)[2] + rnorm(n = 1, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[2])) * bias_sims 
bias_1000_sims <- append(bias_1000_sims, sims)
(sims <- ifelse(is.na(sims), 1, sims))
bias_1000_sims_list[[i]] <- sims
}
(bias_1000_sims <- ifelse(is.na(bias_1000_sims), 1, bias_1000_sims))
bias_1000_sims_list[[10]]
## Load regression HA ~ BWDO
load(file = "./data/DOvsArea_se18.RData")
####
## Create empty dfs to store results
bwdo_out <- data.frame()
ha_out <- data.frame()
years_of_interest <- years_of_calibration
N <- 1000
########################################################################################
## Get simulations with parameter and total uncertainty
cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)
foreach::foreach(j = seq_along(years_of_interest),
                 .packages = c('dplyr', 'lubridate', 'tidyselect', 'tibble'),
                 .export = c('p_funs95', 'get_se_fit_e', 'get_se_fit_w', 
                             'get_bwdo_sum', 'get_ha_sum', 'get_bwdo_ha_par_unc',
                             'one_year_run_vector')) %dopar% {
  (yy <- years_of_interest[j])
  hydro_met_list <- readRDS(file = paste0("./data//hydrometeo_sims_01-21_paral/", yy, ".rds"))
  out_par <- data.frame()
  out_au <- data.frame()
  
  for (x in seq_along(1:N)) {
    
    parvar <- par_all[x, ]
    # Run model with one parameter set, but with all hydrometeorology
    mod_run <- do.call(dplyr::bind_rows, lapply(hydro_met_list, function(y) one_year_run_vector(data = y, params = parvar)))
    
    # Extract only parameter uncertainty
    out_par_pre <- mod_run %>%
      dplyr::group_by(Date) %>%
      dplyr::summarise(bwdoW_par_unc = mean(bwdoW, na.rm = T),
                       bwdoE_par_unc = mean(bwdoE, na.rm = T)) 
      out_par <- dplyr::bind_rows(out_par, out_par_pre)
    # Save all simulations
      out_au <- dplyr::bind_rows(out_au, mod_run)
      print(x)
  }
  
  bias_sims <- seq(1, 122)
  (bias_sims <- ifelse(bias_sims > 30, NA, bias_sims))
  (sim_num_par <- nrow(out_par) / 122)
  (sim_num_all <- nrow(out_au) / 122)
  ###
  bias_1000_sims_all <- rep(1, sim_num_all)  + 
    (rep(coef(mod_bias)[1], sim_num_all) + 
       rnorm(n = sim_num_all, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[1])) +
    (rep(coef(mod_bias)[2], sim_num_all) + 
       rnorm(n = sim_num_all, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[2])) * rep(bias_sims, sim_num_all) 
  (bias_1000_sims_all <- ifelse(is.na(bias_1000_sims_all), 1, bias_1000_sims_all))
  ###
  bias_1000_sims_par <- rep(1, sim_num_par)  + 
    (rep(coef(mod_bias)[1], sim_num_par) + 
       rnorm(n = sim_num_par, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[1])) +
    (rep(coef(mod_bias)[2], sim_num_par) + 
       rnorm(n = sim_num_par, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[2])) * rep(bias_sims, sim_num_par)
  (bias_1000_sims_par <- ifelse(is.na(bias_1000_sims_par), 1, bias_1000_sims_par))
  
  temp_bwdo_old <- get_bwdo_sum_old(indata = out_au, bias_input = bias_vector)
  temp_bwdo <- get_bwdo_sum(indata = out_au, bias_input = bias_1000_sims_all)
  temp_ha_old   <- get_ha_sum_old(indata = out_au, bias_input = bias_vector) 
  temp_ha   <- get_ha_sum(indata = out_au, bias_input = bias_1000_sims_all) 
  
  par_out_dat_old <- get_bwdo_ha_par_unc_old(indata = out_par, bias_input = bias_vector)
  par_out_dat <- get_bwdo_ha_par_unc(indata = out_par, bias_input = bias_1000_sims_par)
  
  glimpse(par_out_dat)
  out_ha_all <- dplyr::left_join(temp_ha, par_out_dat %>% 
                                          dplyr::select(-dplyr::starts_with("bwdo")),
                          by = "Date")
  out_ha_all_old <- dplyr::left_join(temp_ha_old, par_out_dat_old %>% 
                                   dplyr::select(-dplyr::starts_with("bwdo")),
                                 by = "Date")
  
  saveRDS(out_ha_all, paste0("./data/all_sim_09-09_paral/ha_", yy, ".rds"))
  saveRDS(temp_bwdo,  paste0("./data/all_sim_09-09_paral/bwdo_", yy, ".rds"))
  saveRDS(out_ha_all_old, paste0("./data/all_sim_09-09_paral/ha_old_", yy, ".rds"))
  saveRDS(temp_bwdo_old,  paste0("./data/all_sim_09-09_paral/bwdo_old_", yy, ".rds"))
  print(yy)
  rm(out_au); gc()
}
parallel::stopCluster(cl)

##### This is analysis of the results
glimpse(temp_bwdo)
removeit <- bind_cols(temp_bwdo_old %>% select(Date, old = bwdoW_ba_mean, old_nobias = bwdoW_mean, 
                                               lowold = `bwdoW_ba_2.5%`, highold = `bwdoW_ba_97.5%`),
                      temp_bwdo %>% select(new = bwdoW_ba_mean, low = `bwdoW_ba_2.5%`, high = `bwdoW_ba_97.5%`)) %>% 
  mutate(month = lubridate::month(Date))
ggplot(data = removeit, aes(old, new))+
  geom_ribbon(data = removeit, aes(x=old, ymin = low, ymax = high), color = "grey80", alpha = 0.5)+
  geom_ribbon(data = removeit, aes(x=old, ymin = lowold, ymax = highold), color = "grey30", alpha = 0.5)+
  geom_point()+
  facet_wrap(~ month, scale = "free_x")+
  geom_abline(slope = 1)
#####
glimpse(temp_ha)
removeit <- bind_cols(temp_ha_old %>% select(Date, old = ha_w_ParHMUnc_ba_mean, old_nobias = ha_w_ParHMUnc_mean, 
                                               lowold = `ha_w_ParHMUnc_ba_2.5%`, highold = `ha_w_ParHMUnc_ba_97.5%`),
                      temp_ha %>% select(new = ha_w_ParHMUnc_ba_mean, low = `ha_w_ParHMUnc_ba_2.5%`, high = `ha_w_ParHMUnc_ba_97.5%`)) %>% 
  mutate(month = lubridate::month(Date))
ggplot(data = removeit, aes(old, new))+
  geom_ribbon(data = removeit, aes(x=old, ymin = low, ymax = high), color = "grey80", alpha = 0.5)+
  geom_ribbon(data = removeit, aes(x=old, ymin = lowold, ymax = highold), color = "grey30", alpha = 0.5)+
  geom_point()+
  facet_wrap(~ month, scale = "free_x")+
  geom_abline(slope = 1)
ggplot(data = removeit, aes(Date, new))+
  geom_ribbon(data = removeit, aes(x=Date, ymin = low, ymax = high), fill = "grey80", alpha = 0.5)+
  geom_ribbon(data = removeit, aes(x=Date, ymin = lowold, ymax = highold), fill = "grey30", alpha = 0.5)+
  geom_point()+
  facet_wrap(~ month, scale = "free_x")
ggplot(data = removeit, aes(low, lowold))+
  geom_point()
glimpse(temp_ha)
mean((temp_ha$ha_w_AllUnc_var[1:30] - temp_ha_old$ha_w_AllUnc_var[1:30]) / temp_ha_old$ha_w_AllUnc_var[1:30] * 100)
mean(par_out_dat_old$ha_w_Par_ba_var / par_out_dat$ha_w_Par_ba_var)
ggplot()+
  geom_ribbon(data = par_out_dat_old %>% 
                mutate(month = lubridate::month(Date)), 
              aes(x=Date, ymin = `ha_w_Par_ba_2.5%`, ymax = `ha_w_Par_ba_97.5%`), fill = "grey80", alpha = 0.5)+
  geom_ribbon(data = par_out_dat %>% 
                mutate(month = lubridate::month(Date)), 
              aes(x=Date, ymin = `ha_w_Par_ba_2.5%`, ymax = `ha_w_Par_ba_97.5%`), fill = "grey30", alpha = 0.5)+
  geom_point()
ggplot()+
  geom_ribbon(data = temp_ha_old %>% 
                mutate(month = lubridate::month(Date)), 
              aes(x=Date, ymin = `ha_w_AllUnc_ba_2.5%`, ymax = `ha_w_AllUnc_ba_97.5%`), color = "grey80", alpha = 0.5)+
  geom_ribbon(data = temp_ha %>% 
                mutate(month = lubridate::month(Date)), 
              aes(x=Date, ymin = `ha_w_AllUnc_ba_2.5%`, ymax = `ha_w_AllUnc_ba_97.5%`), color = "grey30", alpha = 0.5)+
  geom_point()
#######################################################################################################
## Nutrient Loading and Hydrometeo uncertainty 
cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)
foreach::foreach(j = seq_along(years_of_interest),
                 .packages = c('dplyr', 'lubridate', 'tidyselect', 'tibble'),
                 .export = c('p_funs95', 'get_se_fit_e', 'get_se_fit_w', 
                             'get_bwdo_sum', 'get_ha_sum', 'get_bwdo_ha_par_unc',
                             'one_year_run_vector')) %dopar% {
yy <- years_of_interest[j]
hydro_met_list <- readRDS(file = paste0("./data//hydrometeo_sims_01-21_paral/", yy, ".rds"))
out_hm <- data.frame()
for (k in seq_along(hydro_met_list)) {
  tempdata <- hydro_met_list[[k]]
  tempsims <- do.call(dplyr::bind_rows, apply(par_all, 1, function(x) one_year_run_vector(data = tempdata, params = x)))
  out_hm_pre <- tempsims %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(bwdoW_hm_unc = mean(bwdoW, na.rm = T),
                     bwdoE_hm_unc = mean(bwdoE, na.rm = T)) 
  out_hm <- dplyr::bind_rows(out_hm, out_hm_pre)
}

bias_sims <- seq(1, 122)
(bias_sims <- ifelse(bias_sims > 30, NA, bias_sims))
(sim_num_hm <- nrow(out_hm) / 122)
###
bias_1000_sims_hm <- rep(1, sim_num_hm)  + 
  (rep(coef(mod_bias)[1], sim_num_hm) + 
     rnorm(n = sim_num_hm, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[1])) +
  (rep(coef(mod_bias)[2], sim_num_hm) + 
     rnorm(n = sim_num_hm, mean = 0, sd = sqrt(diag(vcov(mod_bias)))[2])) * rep(bias_sims, sim_num_hm) 
(bias_1000_sims_hm <- ifelse(is.na(bias_1000_sims_hm), 1, bias_1000_sims_hm))
###

hm_out_dat_old <- get_bwdo_ha_hm_unc_old(indata = out_hm, bias_input = bias_vector)
hm_out_dat <- get_bwdo_ha_hm_unc(indata = out_hm, bias_input = bias_1000_sims_hm)
saveRDS(hm_out_dat, paste0("./data/all_sim_09-09_paral/hm_", yy, ".rds"))
saveRDS(hm_out_dat_old, paste0("./data/all_sim_09-09_paral/hm_old_", yy, ".rds"))
}
parallel::stopCluster(cl)
#######################################################################################################
## hindcast error 
res <- data.frame()
for (i in seq_along(years_of_calibration)) {
  tem <- do.call(dplyr::bind_rows, apply(par_all, 1, function(x) one_year_run_vector(data = InpDaif %>% dplyr::filter(year == years_of_calibration[i]), params = x)))
  res <- bind_rows(res, tem)
}
dim(res)
glimpse(res)
rest <- res
rest_bwdo_old <- get_bwdo_sum_old(indata = rest, bias_input = bias_vector) 
glimpse(rest_bwdo_old)
(save_hind_bwdo <- rest_bwdo_old %>% 
  dplyr::select(Date, bwdoE_hind = bwdoE_e_mean, bwdoW_hind = bwdoW_e_ba_mean))
### Write a CSV with hindcasted BWDO
write_csv(save_hind_bwdo, "./data/hindcast_bwdo.csv")
#######################################################################################################
bias_sims <- seq(1, 122)
(bias_sims <- ifelse(bias_sims > 30, NA, bias_sims))
(sim_num <- rep(bias_sims, (nrow(rest) / 122 / 32)))
length(sim_num)
bias_sims_all <- rep(1, length(sim_num))  + 
  (rep(coef(mod_bias)[1], length(sim_num)) + 
     rnorm(n = length(sim_num), mean = 0, sd = sqrt(diag(vcov(mod_bias)))[1])) +
  (rep(coef(mod_bias)[2], length(sim_num)) + 
     rnorm(n = length(sim_num), mean = 0, sd = sqrt(diag(vcov(mod_bias)))[2])) * rep(bias_sims, length(sim_num) / 122) 
(bias_sims_all <- ifelse(is.na(bias_sims_all), 1, bias_sims_all))
(bias_sims_all <- as.numeric(rep(bias_sims_all, 32)))
length(bias_sims_all)
### Run it
rest_ha_old <- get_ha_sum_old(indata = rest, bias_input = bias_vector) 
glimpse(rest_ha_old)
rest_ha <- get_ha_sum(indata = rest, bias_input = bias_sims_all) 
glimpse(rest_ha)
#####################################################################################
rest_ha_old %>% 
  mutate(monthday = format(as.Date(Date), "%m-%d")) %>%
  dplyr::select(Date, monthday, ha_ParHM_ba_var, ha_ParHMMOD_ba_var) %>% 
  dplyr::bind_cols(., rest_ha %>% dplyr::select(ha_AllUnc_ba_var)) %>%
  group_by(monthday) %>% 
  summarise(`Parameter` = mean(ha_ParHM_ba_var),
            `Parameter and Residual error` = mean(ha_ParHMMOD_ba_var),
            `Parameter, Residual error, Transformation` = mean(ha_AllUnc_ba_var)) %>% 
  melt(., id.vars = "monthday") %>% 
  ggplot(., aes(monthday, value, color = variable, group = variable))+
  geom_line()+
  theme(legend.position=c(.80,.87))+
  labs(x = "", y = expression((km)^{4}), color = "Variance")+
  scale_x_discrete(breaks = c("06-01", "06-15", 
                              "07-01", "07-15", 
                              "08-01", "08-15", 
                              "09-01", "09-15"))
########################################################################################
rest_ha %>% 
  mutate(monthday = format(as.Date(Date), "%m-%d")) %>%
  group_by(monthday) %>% 
  summarise(`Parameter` = mean(ha_ParHM_ba_var),
            `Parameter and Residual error` = mean(ha_ParHMMOD_ba_var),
            `Parameter, Residual error, Transformation` = mean(ha_AllUnc_ba_var)) %>% 
  melt(., id.vars = "monthday") %>% 
  ggplot(., aes(monthday, value, color = variable, group = variable))+
  geom_line()+
  theme(legend.position=c(.80,.87))+
  labs(x = "", y = expression((km)^{4}), color = "Variance")+
  scale_x_discrete(breaks = c("06-01", "06-15", 
                              "07-01", "07-15", 
                              "08-01", "08-15", 
                              "09-01", "09-15"))
### Save mean HA hindcast, total and for 2 sections
(save_hind_ha <- rest_ha %>% 
  dplyr::select(Date, ha_e_AllUnc_mean, ha_w_AllUnc_ba_mean))
### Write a CSV
write_csv(save_hind_ha, "./data/hindcast_ha.csv")
### Save variance
write_csv(rest_ha %>% dplyr::select(Date, `Parameter, Residual error, Transformation` = ha_AllUnc_ba_var), "./data/hindcast_var.csv")
#######################################################################################################