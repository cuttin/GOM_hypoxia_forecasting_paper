enddate <- function() "2018-09-30"
startdate <- function() "1997-01-06"
## Get simulations from lm object and new X
simulateX <- function(object, nsim = 1, seed = NULL, X, ...) {
  object$fitted.values <- X
  simulate(object = object, nsim = nsim, seed = seed, ...)
}
## R-squared
E <- function(pred, obs)
{ 
  data <- bind_cols(pred = pred, obs = obs)
  data <- na.omit(data)
  SSres <- sum((data$obs - data$pred) ^ 2)
  SStot <- sum((data$obs - mean(data$obs)) ^ 2)  
  return(1 - SSres/SStot)
}
##
rmse <- function(pred, obs) {
  data <- bind_cols(pred = pred, obs = obs)
  data <- na.omit(data)
  sqrt(mean((data$pred - data$obs) ^ 2))
  }
##
cuberoot <- function(x) sign(x) * abs(x) ^ (1 / 3)
##
MyMerge <- function(x, y) {
  df <- merge(x, y, by = "Date", all.x = TRUE, all.y = TRUE)
  return(df)
}

# Define 4 seasons
define_four_seasons <- function(data) {
  data$season <- NA
  data <- mutate(data, month = month(Date))
  data$season <- ifelse(data$month %in% c(1, 2, 12), "Winter", 
                        ifelse(data$month %in% c(3, 4, 5), "Spring",
                               ifelse(data$month %in% c(6, 7, 8), "Summer", "Fall")))
  data$season <- factor(data$season, levels = c("Spring", "Summer", "Fall", "Winter"))
  data$month <- NULL
  return(data)
}

# Define 4 seasons from months
define_four_seasons_from_month <- function(data) {
  data$season <- NA
  #data <- mutate(data, month = month(Date))
  data$season <- ifelse(data$month %in% c(1, 2, 12), "Winter", 
                        ifelse(data$month %in% c(3, 4, 5), "Spring",
                               ifelse(data$month %in% c(6, 7, 8), "Summer", "Fall")))
  data$season <- factor(data$season, levels = c("Spring", "Summer", "Fall", "Winter"))
  data$month <- NULL
  return(data)
}

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

### vector one year run
one_year_run_vector <- function(data, params) {
  data <- data %>% tibble::add_column(bwdoE = NA, bwdoW = NA, bwdoW_e = NA, bwdoE_e = NA)
  
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
  
  data$bwdoW_e <- data$bwdoW + rnorm(nrow(data), 0, params["sy.W"])
  data$bwdoE_e <- data$bwdoE + rnorm(nrow(data), 0, params["sy.E"])
  data <- data %>% dplyr::select(Date, bwdoW, bwdoW_e, bwdoE, bwdoE_e)
  return(data)
}

###################################################################################################
###################################################################################################

## Get error from the regressions
get_se_fit_e <- function(dat = dat, md = lmHO.E) {
  load(file = "./data/DOvsArea_se18.RData")
  md$fitted.values <- predict(md, data.frame(muDO.E = dat))
  return(unname(predict(md, data.frame(muDO.E = dat), se.fit = TRUE)$se.fit))
}
get_se_fit_w <- function(dat = dat, md = lmHO.W) {
  load(file = "./data/DOvsArea_se18.RData")
  md$fitted.values <- predict(md, data.frame(muDO.W = dat))
  return(unname(predict(md, data.frame(muDO.W = dat), se.fit = TRUE)$se.fit))
} 

###############################################################################
## get bwdo, ha means and quantiles  
## No uncertainty in bias correction
get_bwdo_sum_old <- function(indata, bias_input) {
  indata$bwdoW_ba   <- indata$bwdoW * rep(bias_input, times = (nrow(indata) / 122))
  indata$bwdoW_e_ba <- indata$bwdoW_e * rep(bias_input, times = (nrow(indata) / 122))
  outall_m <- indata %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise_all(funs(mean(., na.rm = TRUE), var(., na.rm = TRUE), !!!p_funs95))
  return(outall_m)
}

## get bwdo, ha means and quantiles 
get_bwdo_sum <- function(indata, bias_input) {
  indata$bwdoW_ba   <- indata$bwdoW * bias_input
  indata$bwdoW_e_ba <- indata$bwdoW_e * bias_input
  outall_m <- indata %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise_all(funs(mean(., na.rm = TRUE), var(., na.rm = TRUE), !!!p_funs95))
  return(outall_m)
}

###################################################################################################
###################################################################################################
## Get Parameter uncertainty
## no uncertainty in bias
get_bwdo_ha_par_unc_old <- function(indata, bias_input) {
  indata$bwdoW_par_unc_ba <- indata$bwdoW_par_unc * rep(bias_input, times = (nrow(indata) / 122))
  
  ha_e_Par <- lmHO.E$coef[1] + lmHO.E$coef[2] * indata$bwdoE_par_unc + lmHO.E$coef[3] * indata$bwdoE_par_unc ^ 2
  ha_w_Par <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_par_unc + lmHO.W$coef[3] * indata$bwdoW_par_unc ^ 2
  ha_w_Par_ba <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_par_unc_ba + lmHO.W$coef[3] * indata$bwdoW_par_unc_ba ^ 2
  
  ha_Par = ha_e_Par + ha_w_Par
  ha_Par_ba = ha_e_Par + ha_w_Par_ba
  
  step_before_outdata <- dplyr::bind_cols(Date = indata$Date,
                                          bwdoW_par_unc    = indata$bwdoW_par_unc,
                                          bwdoE_par_unc    = indata$bwdoE_par_unc,
                                          bwdoW_par_unc_ba = indata$bwdoW_par_unc_ba,
                                          ha_e_Par         = ha_e_Par,
                                          ha_w_Par         = ha_w_Par,
                                          ha_w_Par_ba      = ha_w_Par_ba,  
                                          ha_Par           = ha_e_Par + ha_w_Par,
                                          ha_Par_ba        = ha_e_Par + ha_w_Par_ba)
  out <- step_before_outdata %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise_all(funs(mean(., na.rm = TRUE), var(., na.rm = TRUE), !!!p_funs95))
}

## Get Parameter uncertainty
get_bwdo_ha_par_unc <- function(indata, bias_input) {
  indata$bwdoW_par_unc_ba <- indata$bwdoW_par_unc * bias_input
  
  ha_e_Par <- lmHO.E$coef[1] + lmHO.E$coef[2] * indata$bwdoE_par_unc + lmHO.E$coef[3] * indata$bwdoE_par_unc ^ 2
  ha_w_Par <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_par_unc + lmHO.W$coef[3] * indata$bwdoW_par_unc ^ 2
  ha_w_Par_ba <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_par_unc_ba + lmHO.W$coef[3] * indata$bwdoW_par_unc_ba ^ 2
  
  ha_Par = ha_e_Par + ha_w_Par
  ha_Par_ba = ha_e_Par + ha_w_Par_ba
  
  step_before_outdata <- dplyr::bind_cols(Date = indata$Date,
                                          bwdoW_par_unc    = indata$bwdoW_par_unc,
                                          bwdoE_par_unc    = indata$bwdoE_par_unc,
                                          bwdoW_par_unc_ba = indata$bwdoW_par_unc_ba,
                                          ha_e_Par         = ha_e_Par,
                                          ha_w_Par         = ha_w_Par,
                                          ha_w_Par_ba      = ha_w_Par_ba,  
                                          ha_Par           = ha_e_Par + ha_w_Par,
                                          ha_Par_ba        = ha_e_Par + ha_w_Par_ba)
  out <- step_before_outdata %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise_all(funs(mean(., na.rm = TRUE), var(., na.rm = TRUE), !!!p_funs95))
}

########################################################################################
## Get input uncertainty
## no bias uncertainty
get_bwdo_ha_hm_unc_old <- function(indata, bias_input) {
  indata$bwdoW_hm_unc_ba <- indata$bwdoW_hm_unc * rep(bias_input, times = (nrow(indata) / 122))
  
  ha_e_hm <- lmHO.E$coef[1] + lmHO.E$coef[2] * indata$bwdoE_hm_unc + lmHO.E$coef[3] * indata$bwdoE_hm_unc ^ 2
  ha_w_hm <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_hm_unc + lmHO.W$coef[3] * indata$bwdoW_hm_unc ^ 2
  ha_w_hm_ba <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_hm_unc_ba + lmHO.W$coef[3] * indata$bwdoW_hm_unc_ba ^ 2
  
  ha_hm = ha_e_hm + ha_w_hm
  ha_hm_ba = ha_e_hm + ha_w_hm_ba
  
  step_before_outdata <- dplyr::bind_cols(Date = indata$Date,
                                          bwdoW_hm_unc    = indata$bwdoW_hm_unc,
                                          bwdoE_hm_unc    = indata$bwdoE_hm_unc,
                                          bwdoW_hm_unc_ba = indata$bwdoW_hm_unc_ba,
                                          ha_e_hm         = ha_e_hm,
                                          ha_w_hm         = ha_w_hm,
                                          ha_w_hm_ba      = ha_w_hm_ba,  
                                          ha_hm           = ha_e_hm + ha_w_hm,
                                          ha_hm_ba        = ha_e_hm + ha_w_hm_ba)
  out <- step_before_outdata %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise_all(funs(mean(., na.rm = TRUE), var(., na.rm = TRUE), !!!p_funs95))
}

## Get input uncertainty
get_bwdo_ha_hm_unc <- function(indata, bias_input) {
  indata$bwdoW_hm_unc_ba <- indata$bwdoW_hm_unc * bias_input
  
  ha_e_hm <- lmHO.E$coef[1] + lmHO.E$coef[2] * indata$bwdoE_hm_unc + lmHO.E$coef[3] * indata$bwdoE_hm_unc ^ 2
  ha_w_hm <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_hm_unc + lmHO.W$coef[3] * indata$bwdoW_hm_unc ^ 2
  ha_w_hm_ba <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_hm_unc_ba + lmHO.W$coef[3] * indata$bwdoW_hm_unc_ba ^ 2
  
  ha_hm = ha_e_hm + ha_w_hm
  ha_hm_ba = ha_e_hm + ha_w_hm_ba
  
  step_before_outdata <- dplyr::bind_cols(Date = indata$Date,
                                          bwdoW_hm_unc    = indata$bwdoW_hm_unc,
                                          bwdoE_hm_unc    = indata$bwdoE_hm_unc,
                                          bwdoW_hm_unc_ba = indata$bwdoW_hm_unc_ba,
                                          ha_e_hm         = ha_e_hm,
                                          ha_w_hm         = ha_w_hm,
                                          ha_w_hm_ba      = ha_w_hm_ba,  
                                          ha_hm           = ha_e_hm + ha_w_hm,
                                          ha_hm_ba        = ha_e_hm + ha_w_hm_ba)
  out <- step_before_outdata %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise_all(funs(mean(., na.rm = TRUE), var(., na.rm = TRUE), !!!p_funs95))
}
########################################################################################
## Get all HA uncertainty
## No variation in bias
get_ha_sum_old <- function(indata, bias_input) { #1000000
  indata$bwdoW_ba   <- indata$bwdoW * rep(bias_input, times = (nrow(indata) / 122))
  indata$bwdoW_e_ba <- indata$bwdoW_e * rep(bias_input, times = (nrow(indata) / 122))
  ##
  ha_e_ParHMUnc <- lmHO.E$coef[1] + lmHO.E$coef[2] * indata$bwdoE + lmHO.E$coef[3] * indata$bwdoE ^ 2
  ha_w_ParHMUnc <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW + lmHO.W$coef[3] * indata$bwdoW ^ 2
  ha_w_ParHMUnc_ba <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_ba + lmHO.W$coef[3] * indata$bwdoW_ba ^ 2
  ## all uncertainty BWDO 
  ha_e_ParHMMODUnc <- lmHO.E$coef[1] + lmHO.E$coef[2] * indata$bwdoE_e + lmHO.E$coef[3] * indata$bwdoE_e ^ 2 
  ha_w_ParHMMODUnc <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_e + lmHO.W$coef[3] * indata$bwdoW_e ^ 2
  ha_w_ParHMMODUnc_ba <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_e_ba + lmHO.W$coef[3] * indata$bwdoW_e_ba ^ 2
  ###
  ## All regression unc
  sefit.E <- get_se_fit_e(md = lmHO.E, dat = indata$bwdoE_e) 
  ha_e_AllUnc <- ha_e_ParHMMODUnc + rnorm(length(ha_e_ParHMMODUnc), mean = 0, sd = sigma(lmHO.E)) +
    rnorm(length(ha_e_ParHMMODUnc), mean = rep(0, length(sefit.E)), sd = sefit.E)
  rm(sefit.E)
  sefit.W <- get_se_fit_w(md = lmHO.W, dat = indata$bwdoW_e)
  ha_w_AllUnc <- ha_w_ParHMMODUnc + rnorm(length(ha_w_ParHMMODUnc), mean = 0, sd = sigma(lmHO.W)) + 
    rnorm(length(ha_w_ParHMMODUnc), mean = rep(0, length(sefit.W)), sd = sefit.W)
  rm(sefit.W)
  sefit.W.ba <- get_se_fit_w(md = lmHO.W, dat = indata$bwdoW_e_ba)
  ha_w_AllUnc_ba <- ha_w_ParHMMODUnc_ba + rnorm(length(ha_w_ParHMMODUnc_ba), mean = 0, sd = sigma(lmHO.W)) + 
    rnorm(length(ha_w_ParHMMODUnc_ba), mean = rep(0, length(sefit.W.ba)), sd = sefit.W.ba)
  rm(sefit.W.ba)
  ##
  step_before_outdata <- dplyr::bind_cols(Date = indata$Date,
                                   ha_e_ParHMUnc    = ha_e_ParHMUnc,
                                   ha_w_ParHMUnc    = ha_w_ParHMUnc,
                                   ha_w_ParHMUnc_ba = ha_w_ParHMUnc_ba,  
                                   ha_ParHM         = ha_e_ParHMUnc + ha_w_ParHMUnc,
                                   ha_ParHM_ba      = ha_e_ParHMUnc + ha_w_ParHMUnc_ba,
                                   ##
                                   ha_e_ParHMMODUnc    = ha_e_ParHMMODUnc,
                                   ha_w_ParHMMODUnc    = ha_w_ParHMMODUnc,
                                   ha_w_ParHMMODUnc_ba = ha_w_ParHMMODUnc_ba,
                                   ha_ParHMMOD         = ha_e_ParHMMODUnc + ha_w_ParHMMODUnc,
                                   ha_ParHMMOD_ba      = ha_e_ParHMMODUnc + ha_w_ParHMMODUnc_ba,
                                   ##
                                   ha_e_AllUnc    = ha_e_AllUnc, 
                                   ha_w_AllUnc    = ha_w_AllUnc,
                                   ha_w_AllUnc_ba = ha_w_AllUnc_ba,
                                   ha_AllUnc      = ha_e_AllUnc + ha_w_AllUnc,
                                   ha_AllUnc_ba   = ha_e_AllUnc + ha_w_AllUnc_ba)
  out <- step_before_outdata %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise_all(funs(mean(., na.rm = TRUE), var(., na.rm = TRUE), !!!p_funs95))
  return(out)
}

## Get all HA uncertainty
get_ha_sum <- function(indata, bias_input) { #1000000
  indata$bwdoW_ba   <- indata$bwdoW * bias_input
  indata$bwdoW_e_ba <- indata$bwdoW_e * bias_input
  ##
  ha_e_ParHMUnc <- lmHO.E$coef[1] + lmHO.E$coef[2] * indata$bwdoE + lmHO.E$coef[3] * indata$bwdoE ^ 2
  ha_w_ParHMUnc <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW + lmHO.W$coef[3] * indata$bwdoW ^ 2
  ha_w_ParHMUnc_ba <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_ba + lmHO.W$coef[3] * indata$bwdoW_ba ^ 2
  ## all uncertainty BWDO 
  ha_e_ParHMMODUnc <- lmHO.E$coef[1] + lmHO.E$coef[2] * indata$bwdoE_e + lmHO.E$coef[3] * indata$bwdoE_e ^ 2 
  ha_w_ParHMMODUnc <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_e + lmHO.W$coef[3] * indata$bwdoW_e ^ 2
  ha_w_ParHMMODUnc_ba <- lmHO.W$coef[1] + lmHO.W$coef[2] * indata$bwdoW_e_ba + lmHO.W$coef[3] * indata$bwdoW_e_ba ^ 2
  ###
  ## All regression unc
  sefit.E <- get_se_fit_e(md = lmHO.E, dat = indata$bwdoE_e) 
  ha_e_AllUnc <- ha_e_ParHMMODUnc + rnorm(length(ha_e_ParHMMODUnc), mean = 0, sd = sigma(lmHO.E)) +
    rnorm(length(ha_e_ParHMMODUnc), mean = rep(0, length(sefit.E)), sd = sefit.E)
  rm(sefit.E)
  sefit.W <- get_se_fit_w(md = lmHO.W, dat = indata$bwdoW_e)
  ha_w_AllUnc <- ha_w_ParHMMODUnc + rnorm(length(ha_w_ParHMMODUnc), mean = 0, sd = sigma(lmHO.W)) + 
    rnorm(length(ha_w_ParHMMODUnc), mean = rep(0, length(sefit.W)), sd = sefit.W)
  rm(sefit.W)
  sefit.W.ba <- get_se_fit_w(md = lmHO.W, dat = indata$bwdoW_e_ba)
  ha_w_AllUnc_ba <- ha_w_ParHMMODUnc_ba + rnorm(length(ha_w_ParHMMODUnc_ba), mean = 0, sd = sigma(lmHO.W)) + 
    rnorm(length(ha_w_ParHMMODUnc_ba), mean = rep(0, length(sefit.W.ba)), sd = sefit.W.ba)
  rm(sefit.W.ba)
  ##
  step_before_outdata <- dplyr::bind_cols(Date = indata$Date,
                                          ha_e_ParHMUnc    = ha_e_ParHMUnc,
                                          ha_w_ParHMUnc    = ha_w_ParHMUnc,
                                          ha_w_ParHMUnc_ba = ha_w_ParHMUnc_ba,  
                                          ha_ParHM         = ha_e_ParHMUnc + ha_w_ParHMUnc,
                                          ha_ParHM_ba      = ha_e_ParHMUnc + ha_w_ParHMUnc_ba,
                                          ##
                                          ha_e_ParHMMODUnc    = ha_e_ParHMMODUnc,
                                          ha_w_ParHMMODUnc    = ha_w_ParHMMODUnc,
                                          ha_w_ParHMMODUnc_ba = ha_w_ParHMMODUnc_ba,
                                          ha_ParHMMOD         = ha_e_ParHMMODUnc + ha_w_ParHMMODUnc,
                                          ha_ParHMMOD_ba      = ha_e_ParHMMODUnc + ha_w_ParHMMODUnc_ba,
                                          ##
                                          ha_e_AllUnc    = ha_e_AllUnc, 
                                          ha_w_AllUnc    = ha_w_AllUnc,
                                          ha_w_AllUnc_ba = ha_w_AllUnc_ba,
                                          ha_AllUnc      = ha_e_AllUnc + ha_w_AllUnc,
                                          ha_AllUnc_ba   = ha_e_AllUnc + ha_w_AllUnc_ba)
  out <- step_before_outdata %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise_all(funs(mean(., na.rm = TRUE), var(., na.rm = TRUE), !!!p_funs95))
  return(out)
}

###################################################################################################
###################################################################################################

extract_ymd <- function(x) {
  x <- mutate(x, year = year(Date), month = month(Date), day = day(Date))
}

extract_ym <- function(x) {
  x <- mutate(x, year = year(Date), month = month(Date))
}

safe.max = function(invector) {
  na.pct = sum(is.na(invector))/length(invector)
  if (na.pct == 1) {
    return(NA) }
  else {
    return(max(invector, na.rm = TRUE))
  }
}
safe.min = function(invector) {
  na.pct = sum(is.na(invector))/length(invector)
  if (na.pct == 1) {
    return(NA) }
  else {
    return(min(invector, na.rm = TRUE))
  }
}

get_month_year_start_end <- function(data) {
  data <- mutate(data, year = year(Date), month = month(Date)) %>%
    filter(Date < enddate & Date > startdate)
}
# Log-transformation
sysanal.boxcox <- function(data, lambda1 = 1, lambda2 = 1)
{
  if ( lambda1 == 0 )
  {
    return(ifelse(data > -lambda2, log(data + lambda2), NA))
  }
  else
  {
    return(ifelse(data >= -lambda2, ((data + lambda2) ^ lambda1 - 1) / lambda1, NA))
  }
}

sysanal.boxcox.deriv <- function(data, lambda1 = 1, lambda2 = 1)
{
  return(ifelse(data > -lambda2, (data + lambda2) ^ (lambda1 - 1), NA))
}
#
cbbPalette <- function () c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
##
Palette_segments <- function () c('#e66101', '#5e3c99')
Palette_models_123 <- function () c("#CC79A7", "#009E73", "#0072B2")
Linetypes_models_123 <- function () c("dotdash", "dashed", "solid") 
Palette_months <- function () c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")
Palette_variances <- function () c('#e66101', '#fdb863', '#5e3c99')
Linetypes_hyear_123 <- function () c("dotdash", "dashed", "solid")
Palette_hyear_123 <- function () c("#F0E442", "#009E73", "#0072B2")
Linetypes_seasons <- function () c("dashed", "solid", "dotdash", "longdash")
#Palette_seasons <- function () c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")
#Palette_seasons <- function () c('#00429d', '#73a2c6', '#f4777f', '#93003a')
Palette_seasons <- function () c('#00429d', '#4771b2', '#73a2c6', '#a5d5d8')
Shapes_segments <- function () c(5, 16, 17)
Palette_par_distr <- function () c("#009E73", "#CC79A7")
Linetypes_par_distr <- function () c("twodash", "solid") 
## Extract water year
extract_water_yr <- function(dates, start_month = 10) {
  # Convert dates into POSIXlt
  dates.posix = as.POSIXlt(dates)
  # Year offset
  offset = ifelse(dates.posix$mon >= start_month - 1, 1, 0)
  # Water year
  adj.year = dates.posix$year + 1900 + offset
  # Return the water year
  return(adj.year)
}

shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}
