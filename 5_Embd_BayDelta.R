#Energy used by adults migrating through the bay and delta (EMbd), and timing of freshwater entry
#Bay/delta = Golden Gate to Rio Vista (RKM 0), where the river temperature and velocity records start
#Two swimming behaviours: tidal surfing (maintenance metabolism only) or active swimming at 1 BL/s
#Also summarises the date tracked fish passed Rio Vista, used to weight entry timing (script 10)
#Results, Table 2, Supplement 7
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#####packages required
#####
library(MASS) #fitdistr
#####

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
if(!file.exists("Output/FL_Mass_coefficients.csv")) stop("Run 2_FL_Mass.R first")

#read telemetry summary for bay and delta#####
#one row per tagged fish that later reached Freeport, 2022 (n=11) and 2023 (n=56)
#GG_date: first detection past the Golden Gate; RV_date: first detection near, but not past, Rio Vista
#GG_dist/RV_dist: distance from the Golden Gate (km) at those detections
#Temp_mean: mean temperature experienced between the two detections, C
RV <- read.csv("Data/Telemetry_BayDelta_2022_2023.csv")
RV$GG_date <- as.POSIXct(RV$GG_date_UTC, tz="UTC")
RV$RV_date <- as.POSIXct(RV$RV_date_UTC, tz="UTC")
RV$mig_dur_deltaloc <- as.numeric(difftime(RV$RV_date, RV$GG_date, units = "days"))
RV <- subset(RV, mig_dur_deltaloc > 0) #only fish that were detected later
RV$mig_dist_RV <- RV$RV_dist_km - RV$GG_dist_km
RV <- subset(RV, mig_dist_RV > 0) #only fish that actually made it past GG
RV$mig_rate_km_d <- RV$mig_dist_RV / RV$mig_dur_deltaloc

nrow(RV) #67
table(RV$Year) #11, 56
summary(RV$mig_dur_deltaloc) #days to reach Rio Vista, median 17.6
summary(RV$mig_rate_km_d) #km/d
summary(RV$Temp_mean) #C
#####

#date passed Rio Vista: entry into the model domain#####
#used to weight simulated fish by observed entry timing (script 10)
#pooled across years so that differences between years reflect river conditions, not run timing
RV$RV_J <- as.numeric(format(RV$RV_date, "%j"))
summary(subset(RV, Year == 2022)$RV_J)
summary(subset(RV, Year == 2023)$RV_J)
(Entry_mean <- mean(RV$RV_J)) #259.5, Sep 16
(Entry_sd   <- sd(RV$RV_J))   #16.6 d
(Entry_min  <- min(RV$RV_J))  #234, Aug 22, earliest observed entry
#####

#fish mass from fork length (script 2)#####
FL_Mass_coef <- read.csv("Output/FL_Mass_coefficients.csv")
RV$Mass_kg <- exp(FL_Mass_coef$slope*log(RV$Fork_length_cm*10) + FL_Mass_coef$intercept) #FL in mm
summary(RV$Mass_kg) #kg

#body length from fork length, Conrad and Gutmann 1996 #cm
RV$BL_cm <- ifelse(RV$Fork_length_cm < 68, 1.023 + 1.045*RV$Fork_length_cm, 1.488 + 1.032*RV$Fork_length_cm)
RV$Us_1BL <- RV$BL_cm/100 #swim speed of 1 BL/s, m/s
summary(RV$Us_1BL)
#####

#bioenergetics parameters and functions#####
#Martin et al. 2015
j     = 1.358442*(10^-5)    #MJ/mgO2             #to convert from mgO2 to MJ
c_0   = 65.51               #mgO2/h/kg           #Intercept for mass-specific maintenance rate
b     = -0.217              #unitless            #Mass exponent for mass-specific maintenance rate
d     = 0.068               #C^-1                #Temperature exponent for maintenance rate
c_1   = 420.6               #mgO2/h/kg/m^-gamma  #Intercept for swimming cost scalar (1 kg)
beta  = -0.369              #unitless            #Mass exponent for mass-specific swim costs
gamma = 2.130               #unitless            #Velocity exponent for mass-specific swim costs
yr    = 1.82                #unitless            #Recovery yield

Pm.func <- function(M, temp){return(c_0*(M^b)*exp(d*temp))} #mgO2/h/kg #maintenance metabolism
Ps.func <- function(Us, M){return(c_1*(M^beta)*Us^gamma)} #mgO2/h/kg #swimming, Us in m/s
Pc_temp.func <- function(temp){return(274.7371 + 49.0283*temp - 1.4244*temp^2)} #critical power, FitzGerald & Martin 2022
COST.func <- function(Pm, Ps, Pc){return(Pm + Ps + Pm*pmax(0, (yr*(Pm+Ps-Pc))/(Pc-Pm)) + pmax(0, (yr-1)*(Pm+Ps-Pc)))} #mgO2/h/kg
#####

#tidal surfing: fish ride flood tides and rest on ebb tides, so use maintenance metabolism only#####
RV$Pm <- Pm.func(M = RV$Mass_kg, temp = RV$Temp_mean)
RV$Em_baydelta <- RV$Pm*24*j*RV$mig_dur_deltaloc #MJ/kg
summary(RV$Em_baydelta) #mean 1.2 MJ/kg
sd(RV$Em_baydelta)      #0.8 MJ/kg
mean(RV$Em_baydelta)/9.5*100 #percent of EI
#####

#active swimming at 1 BL/s for the whole bay/delta transit#####
RV$Ps <- Ps.func(Us = RV$Us_1BL, M = RV$Mass_kg)
RV$Pc <- Pc_temp.func(temp = RV$Temp_mean)
RV$COST <- COST.func(Pm = RV$Pm, Ps = RV$Ps, Pc = RV$Pc)
RV$Em_baydelta_1BL <- RV$COST*24*j*RV$mig_dur_deltaloc #MJ/kg
summary(RV$Em_baydelta_1BL) #mean 2.2 MJ/kg, roughly twice tidal surfing
#####

#distribution of tidal-surfing energy use, for the survival model#####
#maximum likelihood fits; lowest AIC is the best fit
fitN  <- fitdistr(RV$Em_baydelta, "normal")
fitW  <- fitdistr(RV$Em_baydelta, "weibull", lower = c(0.01, 0.01))
fitG  <- fitdistr(RV$Em_baydelta, "gamma", lower = c(0.01, 0.01))
fitLN <- fitdistr(RV$Em_baydelta, "lognormal")
data.frame(distribution = c("normal", "Weibull", "gamma", "lognormal"),
           AIC = c(AIC(fitN), AIC(fitW), AIC(fitG), AIC(fitLN))) #lognormal lowest
fitLN #meanlog 0.028 (SE 0.072), sdlog 0.590 (SE 0.051)

(Embd_meanlog <- fitLN$estimate["meanlog"])
(Embd_sdlog   <- fitLN$estimate["sdlog"])
exp(Embd_meanlog + Embd_sdlog^2/2) #mean of fitted distribution, MJ/kg

#for the sensitivity analysis (script 10): lognormal with the same sdlog, centred on the
#active-swimming mean of 2.2 MJ/kg
(Embd_meanlog_cruise <- log(2.2) - 0.590^2/2) #0.614

Embd_out <- data.frame(n = nrow(RV),
                       Embd_tidalsurf_mean = mean(RV$Em_baydelta), Embd_tidalsurf_sd = sd(RV$Em_baydelta),
                       Embd_1BL_mean = mean(RV$Em_baydelta_1BL), Embd_1BL_sd = sd(RV$Em_baydelta_1BL),
                       meanlog = unname(Embd_meanlog), sdlog = unname(Embd_sdlog),
                       meanlog_cruise = Embd_meanlog_cruise,
                       Entry_mean_J = Entry_mean, Entry_sd_J = Entry_sd, Entry_min_J = Entry_min)
write.csv(Embd_out, "Output/Embd_and_entry_parameters.csv", row.names = FALSE)
#####

###END###
