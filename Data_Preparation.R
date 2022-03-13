library(zoo)
library(ggplot2)

# Initialize data variables
gdp_data = read.csv("Data/GDP_PERC_ITALY.csv")
btp_data = read.csv("Data/BTP_QUARTERLY_ITALY.csv")
bot_data = read.csv("Data/BOT_QUARTERLY_ITALY.csv")
unempl_data = read.csv("Data/UNEMPL_Q_ITALY.csv")
cpi_data = read.csv("Data/CPI_Q_ITALY.csv")

# Creation of dataframes from the various datasets 
gdp_data$quarter = as.yearqtr(gdp_data$TIME, format="%Y-Q%q")
gdp_data$qvar = as.Date(gdp_data$quarter)
gdp_df = data.frame(date = as.Date(gdp_data$qvar, "%Y-%m-%d"), value=gdp_data$Value)

btp_data$quarter = as.yearqtr(btp_data$TIME, format="%Y-Q%q")
btp_data$qvar = as.Date(btp_data$quarter)
btp_df = data.frame(date = as.Date(btp_data$qvar, "%Y-%m-%d"), value=btp_data$Value)

bot_data$quarter = as.yearqtr(bot_data$TIME, format="%Y-Q%q")
bot_data$qvar = as.Date(bot_data$quarter)
bot_df = data.frame(date = as.Date(bot_data$qvar, "%Y-%m-%d"), value=bot_data$Value)

unempl_data$quarter = as.yearqtr(unempl_data$TIME, format="%Y-Q%q")
unempl_data$qvar = as.Date(unempl_data$quarter)
unempl_df = data.frame(date = as.Date(unempl_data$qvar, "%Y-%m-%d"), value=unempl_data$Value)
unempl_df <- unempl_df[order(unempl_df$date),]

cpi_data$qvar = as.Date(cpi_data$QUARTERS, format="%d/%m/%Y")
cpi_df = data.frame(date = cpi_data$qvar, value=cpi_data$INFLATION)

# Variables for HMM and calculation
gdp_rate <- as.numeric(gdp_df[, 2])
Date <- as.numeric(gdp_df[, 1])

y_gdp <- as.numeric(gdp_rate)
model_gdp <- depmix(y_gdp ~ 1, data=data.frame(y_gdp), nstates=3)

U_rate <- as.numeric(unempl_df[, 2])
Date <- as.numeric(unempl_df[, 1])

y_U <- as.numeric(U_rate)
model_U <- depmix(y_U ~ 1, data=data.frame(y_U), nstates=3)

fmodel_gdp <- fit(model_gdp)
fmodel_gdp
states_mu_sigma_gdp <- summary(fmodel_gdp)
trans_prob_gdp <- fmodel_gdp@trDens

fmodel_U <- fit(model_U)
fmodel_U
states_mu_sigma_U <- summary(fmodel_U)
trans_prob_U <- fmodel_U@trDens

df_trans_prob <- data.frame(c("State 1", "State 2", "State 3"), c(round(trans_prob_gdp[1,1,1], 3), round(trans_prob_gdp[1,1,2], 3), round(trans_prob_gdp[1,1,3], 3)), c(0.057, 0.121, 0.103), round(c(trans_prob_gdp[1,2,1], trans_prob_gdp[1,2,2], trans_prob_gdp[1,2,3]), 3), c(0.031, 0.121, 0.072), round(c(trans_prob_gdp[1,3,1], trans_prob_gdp[1,3,2], trans_prob_gdp[1,3,3]), 3), c(0.072, "NA", 0.129), round(c(trans_prob_U[1,1,1], trans_prob_U[1,1,2], trans_prob_U[1,1,3]), 3), c(0.072, "NA", 0.033), round(c(trans_prob_U[1,2,1], trans_prob_U[1,2,2], trans_prob_U[1,2,3]), 3), c("NA", 0.019, 0.033), round(c(trans_prob_U[1,3,1], trans_prob_U[1,3,2], trans_prob_U[1,3,3]), 3), c(0.072, 0.019, 0.046))
df_emis_distr <- data.frame(c("State 1", "State 2", "State 3"), c(states_mu_sigma_gdp[1,1], states_mu_sigma_gdp[2,1], states_mu_sigma_gdp[3,1]), c(0.056, 0.264, 0.165), c(states_mu_sigma_gdp[1,2], states_mu_sigma_gdp[2,2], states_mu_sigma_gdp[3,2]), c(0.035, 0.177, 0.107), c(states_mu_sigma_U[1,1], states_mu_sigma_U[2,1], states_mu_sigma_U[3,1]), c(0.177, 0.120, 0.109), c(St.deviation = states_mu_sigma_U[1,2], states_mu_sigma_U[2,2], states_mu_sigma_U[3,2]), c(0.127, 0.089, 0.087))

MLEse_gdp = standardError(fmodel_gdp)
MLEse_U = standardError(fmodel_U)

estStates_gdp <- posterior(fmodel_gdp)
estStates_U <- posterior(fmodel_U)

HMM_df_gdp <- data.frame(state = 1:3, 
                         mu = states_mu_sigma_gdp[,1], 
                         sigma = states_mu_sigma_gdp[,2])

df_to_plot_gdp <- estStates_gdp %>% as.data.frame() %>%
  left_join(HMM_df_gdp)

HMM_df_U <- data.frame(state = 1:3, 
                       mu = states_mu_sigma_U[,1], 
                       sigma = states_mu_sigma_U[,2])

df_to_plot_U <- estStates_U %>% as.data.frame() %>%
  left_join(HMM_df_U)

# DLM - Random Walk + Noise

# MLE for unknown variances of a random walk plus noise and model specification with Filtering
gdp_ts <- as.ts(gdp_df[3:96,2])
buildrw <- function(param){dlmModPoly(order=1, dV=param[1], dW=param[2], m0=gdp_ts[1])}
outMLE_rw <- dlmMLE(gdp_ts, parm = rep(100, 2), buildrw, lower=c(0.00001, 0), hessian = TRUE)

gdpMod <- buildrw(outMLE_rw$par)
rw <- dlm(m0=0, C0=10, FF=1, V=outMLE_rw$par[1], GG=1, W=outMLE_rw$par[2])

# Computation of asymptotic standard errors 
hessian_rw <- outMLE_rw$hessian
AsymCov_rw = solve(hessian_rw) 
st_errors_rw <- sqrt(diag(AsymCov_rw))

# Smoothing
outSmooth = dlmSmooth(gdp_ts, gdpMod)
rw_filter_results = dlmFilter(gdp_ts, gdpMod)

SmoothingEstimates = window(outSmooth$s,start=start(gdp_ts)[1]-1)
SmoothingEstimates <- dropFirst(SmoothingEstimates)

# Standard Deviations from Smoothing results
n <- length(gdp_ts)
listS <- dlmSvd2var(outSmooth$U.S, outSmooth$D.S)
sqrtS <- sqrt(unlist(listS))

# Computation of the Credible Intervals
lb_estimates <- SmoothingEstimates - qnorm(0.975)*sqrtS[-1]
ub_estimates <- SmoothingEstimates + qnorm(0.975)*sqrtS[-1]

# Data for GDP plot
GDP <- data.frame(quarters=gdp_df[3:96,1], data=as.numeric(gdp_ts), estimates=SmoothingEstimates, lb=lb_estimates, ub=ub_estimates)
GDP <- GDP %>% gather(data, estimates, -quarters, -lb, -ub)

# Static regression model - DLM

# Regressors with and without BOT
regressors_lag_first <- matrix(c(bot_df[1:94,]$value, btp_df[1:94,]$value, cpi_df[1:94,]$value), ncol = 3)
regressors_lag <- matrix(c(btp_df[1:94,]$value, cpi_df[1:94,]$value), ncol = 2)

# DLM model with regressors with BOT
buildDLM_first <- function(param){dlmModReg(regressors_lag_first, addInt = TRUE, dV=param, dW=c(0,0,0,0))}
fit_first = dlmMLE(gdp_df[3:96,]$value, parm=1, buildDLM_first, lower=c(0.00001), hessian = TRUE)

# Linear regression with regressors with BOT
outLM_first <- lm(gdp_df[3:96,]$value ~ regressors_lag_first)
LM_summary_first <- summary(outLM_first)
LM_residuals_first <- sum(outLM_first$residuals^2)/outLM_first$df.residual

# Building the DLM model with static coefficients and with regressors without BOT
buildDLM <- function(param){dlmModReg(regressors_lag, addInt = TRUE, dV=param, dW=c(0,0,0))}
fit = dlmMLE(gdp_df[3:96,]$value, parm=1, buildDLM, lower=c(0.00001), hessian = TRUE)

modStaticDLM = buildDLM(fit$par)
staticDLM_filter_results = dlmFilter(gdp_ts, modStaticDLM)
results_staticDLM <- staticDLM_filter_results$m[length(gdp_ts)+1,]

# Asymptotic standard errors computation
hessian_staticDLM <- fit$hessian
AsymCov_staticDLM = solve(hessian_staticDLM) 
st_errors_staticDLM <- sqrt(diag(AsymCov_staticDLM))

# Linear regression with regressors without BOT
outLM <- lm(gdp_df[3:96,]$value ~ regressors_lag)
LM_summary <- summary(outLM)
LM_residuals <- sum(outLM$residuals^2)/outLM$df.residual

# Dynamic regression model - DLM

# Build the dynamic DLM model
buildDynamicLM = function(param){dlmModReg(regressors_lag, dV=param[1], dW=param[2:4], m0=rep(0, 3))}
outMLE_dyn = dlmMLE(gdp_ts, rep(0, 4), buildDynamicLM, lower=c(0.000001, 0, 0, 0), hessian=TRUE)

modDynamicLM = buildDynamicLM(outMLE_dyn$par)

# Asymptotic standard errors computation
hessian_DynamicDLM <- outMLE_dyn$hessian
AsymCov_DynamicDLM = solve(hessian_DynamicDLM) 
st_errors_DynamicDLM <- sqrt(diag(AsymCov_DynamicDLM))

# Compute the smoothing and forecasting
outSmoothLM = dlmSmooth(gdp_ts, modDynamicLM)
dynDLM_filter_results = dlmFilter(gdp_ts, modDynamicLM)

# Credible intervals Beta BTP

SmoothingEstimates_btp = dropFirst(outSmoothLM$s[,2])
SmoothingEstimates_cpi = dropFirst(outSmoothLM$s[,3])

listS_btp_sqrt <- dlmSvd2var(outSmoothLM$U.S, outSmoothLM$D.S) %>% lapply(function(x) x[2,2]) %>% unlist() %>% sqrt()
listS_cpi_sqrt <- dlmSvd2var(outSmoothLM$U.S, outSmoothLM$D.S) %>% lapply(function(x) x[3,3]) %>% unlist() %>% sqrt()

# Computation of the Credible Intervals for BTP and CPI
lb_estimates_btp <- SmoothingEstimates_btp - qnorm(0.975)*listS_btp_sqrt[-1]
ub_estimates_btp <- SmoothingEstimates_btp + qnorm(0.975)*listS_btp_sqrt[-1]

lb_estimates_cpi <- SmoothingEstimates_cpi - qnorm(0.975)*listS_cpi_sqrt[-1]
ub_estimates_cpi <- SmoothingEstimates_cpi + qnorm(0.975)*listS_cpi_sqrt[-1]

# Data regarding Betas from smoothing results
beta_btp <- data.frame(time = gdp_df[3:96,]$date, value = dropFirst(outSmoothLM$s[,2]), lb=lb_estimates_btp, ub=ub_estimates_btp) 
#beta_btp <- beta_btp %>% gather(value, -time, -lb, -ub)

beta_cpi <- data.frame(time = gdp_df[3:96,]$date, value = dropFirst(outSmoothLM$s[,3]), lb=lb_estimates_cpi, ub=ub_estimates_cpi)
#beta_cpi <- beta_cpi %>% gather(value, -time, -lb, -ub)

# Errors calculations

# We calculate MAE, MSE and MAPE based on forecasting errors

LR_mae_3 = mean(abs(outLM_first$residuals))
LR_mse_3 = mean((outLM_first$residuals)^2)
LR_mape_3 = mean(abs(outLM_first$residuals / gdp_df[3:96,]$value))
LR_errors_summary_3 = c(LR_mae_3, LR_mse_3, LR_mape_3)

LR_mae_2 = mean(abs(outLM$residuals))
LR_mse_2 = mean((outLM$residuals)^2)
LR_mape_2 = mean(abs(outLM$residuals / gdp_df[3:96,]$value))
LR_errors_summary_2 = c(LR_mae_2, LR_mse_2, LR_mape_2)

rw_mae = mean(abs(gdp_ts - rw_filter_results$f))
rw_mse = mean((gdp_ts - rw_filter_results$f)^2)
rw_mape = mean(abs((gdp_ts - rw_filter_results$f) / gdp_ts))
rw_errors_summary = c(rw_mae, rw_mse, rw_mape)

sdlm_mae = mean(abs(gdp_ts - staticDLM_filter_results$f))
sdlm_mse = mean((gdp_ts - staticDLM_filter_results$f)^2)
sdlm_mape = mean(abs((gdp_ts - staticDLM_filter_results$f) / gdp_ts))
sdlm_errors_summary = c(sdlm_mae, sdlm_mse, sdlm_mape)

dyn_mae = mean(abs(gdp_ts - dynDLM_filter_results$f))
dyn_mse = mean((gdp_ts - dynDLM_filter_results$f)^2)
dyn_mape = mean(abs((gdp_ts - dynDLM_filter_results$f) / gdp_ts))
dyn_errors_summary = c(dyn_mae, dyn_mse, dyn_mape)


