#Calibration of the moving fraction, and observed vs predicted migration rates
#The moving fraction is the fraction of each daily time step that a fish migrates; it rests at
#maintenance metabolism for the remainder. Its distribution among individuals (shape) is calibrated on
#fish tracked in 2012; only its mean is fitted to the 2022-2023 fish, by matching the observed median
#migration rate between RKM 20 and 150. Simulated and observed rates are compared with KS tests
#Methods, Results, Table 3; must run script 7 first
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
if(!file.exists("Output/Landscapes/data_ID_Sacramento.csv")) stop("Run 7_EnergyLandscapes.R first")

#observed migration rates#####
#tagged fish detected at both RKM 20 and RKM 150 of the mainstem Sacramento River (RKM 0 = Rio Vista)
#2012 (n=18) to calibrate the shape of the moving fraction; 2022 (n=5) and 2023 (n=46) to fit its mean
obs <- read.csv("Data/Telemetry_MigrationRate_RKM20_150.csv")
obs$RKM20_date  <- as.POSIXct(obs$RKM20_date_UTC,  tz="UTC")
obs$RKM150_date <- as.POSIXct(obs$RKM150_date_UTC, tz="UTC")
obs$entry <- as.integer(format(obs$RKM20_date, "%j")) #day of year the fish passed RKM 20
obs$transit_d <- as.numeric(difftime(obs$RKM150_date, obs$RKM20_date, units = "days"))
obs$rate <- (150 - 20)/obs$transit_d #km/d
table(obs$Year)

o12 <- subset(obs, Year == 2012)
tst <- subset(obs, Year %in% c(2022, 2023))
(obs_median <- median(tst$rate))             #20.0 km/d
(obs_CV <- sd(tst$rate)/mean(tst$rate))      #0.70

#rows of the landscape at RKM 20 and RKM 150
data_ID <- read.csv("Output/Landscapes/data_ID_Sacramento.csv")
rowA <- which.min(abs(data_ID$RKM - 20))
rowB <- which.min(abs(data_ID$RKM - 150))
#####

#Model migration rate for each fish, swimming continuously, under each critical power formulation#####
#each simulated fish starts at RKM 20 on the day its tracked counterpart passed RKM 20
#and steps through the migration-speed landscape exactly as in the migration simulation (script 9)
arms <- c("PcBase", "PcCoastal")
arm_labels <- c("Base-model critical power", "Aerobic scope, coastal fall-run population")
Table3 <- data.frame(Model = "Observed", n_fish = nrow(tst), free_median_km_d = NA, fitted_mean_moving_fraction = NA,
                     median_km_d = obs_median, CV = obs_CV, KS_D = NA, KS_p = NA)

for(a in 1:length(arms)){
  print(arm_labels[a])
  #migration-speed landscapes (m/s) for each year; columns are day of year
  Um_2012 <- readRDS(paste0("Output/Landscapes/Landscape_", arms[a], "_VEL_10pct_2012.rds"))$df_Um
  Um_2022 <- readRDS(paste0("Output/Landscapes/Landscape_", arms[a], "_VEL_10pct_2022.rds"))$df_Um
  Um_2023 <- readRDS(paste0("Output/Landscapes/Landscape_", arms[a], "_VEL_10pct_2023.rds"))$df_Um
  obs$free_rate <- NA
  for(i in 1:nrow(obs)){
    if(obs$Year[i] == 2012){Um_matrix <- Um_2012}
    if(obs$Year[i] == 2022){Um_matrix <- Um_2022}
    if(obs$Year[i] == 2023){Um_matrix <- Um_2023}
    rown <- rowA
    coln <- obs$entry[i]
    steps <- 0
    while(rown < rowB & trunc(coln) <= ncol(Um_matrix)){
      steps <- steps + 1
      if(steps > 20000) break
      Um <- Um_matrix[trunc(rown), trunc(coln)]
      if(!is.finite(Um) | Um <= 0){ #fish cannot advance; waits a day
        coln <- coln + 1
        next
      }
      kmperday <- Um*60*60*24/1000
      mig_dist_km <- ifelse(kmperday >= 2, 2, kmperday) #fish moves a maximum of one 2 km cell
      mig_time_d <- ifelse((1/kmperday*mig_dist_km) >= 1, 1, (1/kmperday*mig_dist_km)) #and a maximum of one day
      rown <- rown + mig_dist_km/2
      coln <- coln + mig_time_d
    }
    if(rown >= rowB){obs$free_rate[i] <- (150 - 20)/(coln - obs$entry[i])} #NA if the fish did not reach RKM 150
  }

  #moving fraction implied for each 2012 fish = observed rate / continuous-swimming rate
  o12 <- subset(obs, Year == 2012)
  duty12 <- o12$rate/o12$free_rate
  duty12 <- duty12[is.finite(duty12) & duty12 > 0]
  shape <- duty12/mean(duty12) #shape of the distribution, mean 1

  #fit the mean of the moving fraction to the 2022-2023 fish, by matching the observed median
  tst <- subset(obs, Year %in% c(2022, 2023))
  FIT <- median(tst$rate)/median(tst$free_rate, na.rm = TRUE)

  #simulated rates: each 2022-2023 fish given a moving fraction drawn from the 2012 shape, 400 times
  #a fish cannot migrate more than 100% of the time
  set.seed(7)
  sim_rate <- c()
  for(r in 1:400){
    sim_rate <- c(sim_rate, tst$free_rate*pmin(1, FIT*sample(shape, nrow(tst), replace = TRUE)))
  }
  sim_rate <- sim_rate[is.finite(sim_rate)]

  #two-sample Kolmogorov-Smirnov test against the observed rates
  #(ties in the simulated rates give a warning that the p-value is approximate; suppressed)
  KS <- suppressWarnings(ks.test(sim_rate, tst$rate))

  print(c(free_median = median(tst$free_rate, na.rm = TRUE), shape_CV = sd(shape)/mean(shape), FIT = FIT,
          median = median(sim_rate), CV = sd(sim_rate)/mean(sim_rate), D = unname(KS$statistic), p = KS$p.value))
  Table3 <- rbind(Table3, data.frame(Model = arm_labels[a], n_fish = nrow(tst),
                                     free_median_km_d = median(tst$free_rate, na.rm = TRUE), fitted_mean_moving_fraction = FIT,
                                     median_km_d = median(sim_rate), CV = sd(sim_rate)/mean(sim_rate),
                                     KS_D = unname(KS$statistic), KS_p = KS$p.value))

  #keep the calibration of the adopted model for the simulations (script 9)
  if(arms[a] == "PcCoastal"){
    shape_PcCoastal <- shape
    FIT_PcCoastal <- FIT
  }
}
#####

#TABLE 3#####
Table3
write.csv(Table3, "Output/Table3_MigrationRate_calibration.csv", row.names = FALSE)
#####

#Moving fraction levels for the simulations#####
#the distribution is represented by nine quantiles of the 2012 shape, scaled by the fitted mean;
#each level is simulated separately and survival is averaged across levels (scripts 9 and 10)
#sensitivity analysis (Table 4): fitted mean -25% and +25%
quantiles <- seq(0.05, 0.95, length.out = 9)
MovingFraction <- data.frame(quantile = quantiles,
                             shape_quantile = unname(quantile(shape_PcCoastal, quantiles)),
                             fitted_mean = FIT_PcCoastal)
MovingFraction$level      <- pmin(1, MovingFraction$shape_quantile*FIT_PcCoastal)
MovingFraction$level_m25  <- pmin(1, MovingFraction$shape_quantile*FIT_PcCoastal*0.75)
MovingFraction$level_p25  <- pmin(1, MovingFraction$shape_quantile*FIT_PcCoastal*1.25)
MovingFraction
write.csv(MovingFraction, "Output/MovingFraction_levels.csv", row.names = FALSE)
saveRDS(MovingFraction, "Output/MovingFraction_levels.rds") #full numerical precision, for script 9
#####

###END###
