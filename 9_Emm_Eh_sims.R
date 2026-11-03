#Energy use for adults during migration (EMr) and holding (EH) on the mainstem Sacramento River
#Simulates one fish entering freshwater (passing Rio Vista, RKM 0) on every day of year 180-365 in 2022
#and 2023, migrating to the mean spawning location and holding until spawning
#Run for each level of the moving fraction (script 8) and for each sensitivity scenario (Table 4)
#Results, Figures 4 and 5, Table 4; must run scripts 7 and 8 first
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
if(!file.exists("Output/MovingFraction_levels.rds")) stop("Run 8_MovingFraction_Calibration.R first")
dir.create("Output/Sims", showWarnings = FALSE, recursive = TRUE)

starttime <- Sys.time()

#read in data#####
data_ID <- read.csv("Output/Landscapes/data_ID_Sacramento.csv") #RKM of each landscape row
MovingFraction <- readRDS("Output/MovingFraction_levels.rds")
MovingFraction

#scenarios: landscape (script 7) and moving fraction levels
#Baseline: aerobic scope of the coastal fall-run population, 10th percentile velocity, fitted moving fraction
#the other scenarios each change one of these for the sensitivity analysis (Table 4)
#Figure: one simulated fish per entry date at the fitted mean moving fraction (Figures 4 and 5)
scenarios <- list(
  Baseline            = list(landscape = "PcCoastal_VEL_10pct",               duty = MovingFraction$level),
  VEL_MEAN            = list(landscape = "PcCoastal_VEL_MEAN",                duty = MovingFraction$level),
  PejusMinus1SD       = list(landscape = "PcCoastal_pejusMinus1SD_VEL_10pct", duty = MovingFraction$level),
  PejusPlus1SD        = list(landscape = "PcCoastal_pejusPlus1SD_VEL_10pct",  duty = MovingFraction$level),
  MovingFraction_m25  = list(landscape = "PcCoastal_VEL_10pct",               duty = MovingFraction$level_m25),
  MovingFraction_p25  = list(landscape = "PcCoastal_VEL_10pct",               duty = MovingFraction$level_p25),
  Figure              = list(landscape = "PcCoastal_VEL_10pct",               duty = MovingFraction$fitted_mean[1]))
#####

#####Salmon distributions#####
#phenology is in Julian day; spatial distribution in m upstream of RKM 0
#model fits from years of redd survey data (FitzGerald et al.)
set.seed(789)
Spawning <- rlnorm(1000, meanlog=3.78613934, sdlog=0.04603214)*7 #spawn date; lognormal on a week scale, mean ~Nov 4
spatialreddsdistrib <- rweibull(1000, shape=21.15206*1.609, scale=292.45399*1.609)*1000 #spawn location
spatialreddsdistrib <- ifelse(spatialreddsdistrib > max(data_ID$m_ds_Keswick), max(data_ID$m_ds_Keswick), spatialreddsdistrib) #so, cannot pick spawn site us of Keswick

SpawnSite_goal <- mean(spatialreddsdistrib) #m upstream of RKM 0, ~457.3 km
goal_row  <- which.min(abs(data_ID$RKM - SpawnSite_goal/1000)) #row of the spawning grounds
start_row <- which.min(abs(data_ID$RKM - 0)) #always start at Rio Vista
#####

#SIMULATION OF ENERGY USE DURING MIGRATION#####
tempres <- 1    #temporal resolution, 1 day
spatialres <- 2 #spatial resolution, 2 km

#steps one fish through the landscape from its entry date until it reaches the spawning grounds
#each step the fish moves at most one 2 km cell and at most one day
#moving fraction (duty): the fish migrates for that fraction of each step and rests at maintenance
#metabolism (Pm) for the remainder, so ground speed is scaled by duty and the cost of the step is
#duty*COST + (1-duty)*Pm, where COST is energy per unit elapsed time while migrating (script 7)
#where aerobic scope has fallen to maintenance the migration speed is 0: the fish cannot advance,
#and holds a full day at maintenance metabolism
EnergyMig.func <- function(entrydate, duty){
  rown <- start_row
  coln <- entrydate
  COST_sum <- 0 #mgO2/kg
  mig_time_cumsum <- 0 #d
  days_halted <- 0
  MigDur <- c()
  RKM <- c()
  temp <- c()
  Uc <- c()
  while(rown < goal_row & trunc(coln) <= ncol(df_t)){
    coln_rd <- trunc(coln) #trunc rounds down to nearest integer
    rown_rd <- trunc(rown)
    RKM    <- c(RKM, rown*spatialres - 2) #so, start at RV is row 1, RKM 0
    MigDur <- c(MigDur, coln_rd - entrydate)
    temp   <- c(temp, df_t[rown_rd, coln_rd]) #C
    Uc     <- c(Uc, df_v[rown_rd, coln_rd])   #m/s
    Um     <- df_Um[rown_rd, coln_rd]         #m/s
    Pm     <- df_Pm[rown_rd, coln_rd]         #mgO2/h/kg
    COST   <- df_COST[rown_rd, coln_rd]       #mgO2/h/kg
    kmperday <- Um*60*60*24/1000*duty #convert Um to km/d, scaled by the moving fraction

    #spatial movement
    mig_dist_km <- ifelse(kmperday >= spatialres, spatialres, kmperday) #fish will move a max of spatialres
    mig_dist_move <- mig_dist_km/spatialres #how many cells fish will move, max of 1
    #temporal movement and energy used
    if(mig_dist_move == 0){ #fish cannot advance; holds a day at maintenance metabolism
      mig_time_d <- tempres
      step_rate <- Pm
      days_halted <- days_halted + 1
    }else{
      mig_time_d <- ifelse((1/kmperday*mig_dist_km) >= tempres, tempres, (1/kmperday*mig_dist_km))
      step_rate <- duty*COST + (1 - duty)*Pm
    }
    mig_time_cumsum <- mig_time_cumsum + mig_time_d
    COST_sum <- COST_sum + step_rate*mig_time_d*24 #mgO2/kg/h x hours in step

    rown <- rown + mig_dist_move
    coln <- coln + mig_time_d
    if(coln >= ncol(df_t)) break #end of temperature record
  }
  return(data.frame(EntryDate_J = entrydate, MigDur_d = max(MigDur), RKM = max(RKM),
                    Temp_mean = mean(temp), Uc_mean = mean(Uc), COST = COST_sum,
                    mig_time_d = mig_time_cumsum, Days_halted = days_halted))
}
#####

#run simulations#####
#simulate at every entry date
EntryDates <- seq(180, 365, 1)
nspawn <- 500 #spawn dates sampled per fish

for(s in 1:length(scenarios)){
  for(lev in 1:length(scenarios[[s]]$duty)){
    duty <- scenarios[[s]]$duty[lev]
    MH <- data.frame()
    for(year in c(2022, 2023)){
      landscape <- readRDS(paste0("Output/Landscapes/Landscape_", scenarios[[s]]$landscape, "_", year, ".rds"))
      df_t    <- landscape$df_t
      df_v    <- landscape$df_v
      df_Um   <- landscape$df_Um
      df_Pm   <- landscape$df_Pm
      df_COST <- landscape$df_COST

      for(i in 1:length(EntryDates)){
        fish <- EnergyMig.func(entrydate = EntryDates[i], duty = duty)

        #did fish reach spawning grounds
        fish$Arrival_J <- fish$EntryDate_J + fish$mig_time_d #when did fish arrive to spawning grounds
        fish$SpawnSitediff <- fish$RKM - SpawnSite_goal/1000 #how close was fish to spawning grounds
        fish$SpawnSiteSuccess <- 1 #Default: assume that fish reached spawning grounds
        if(fish$MigDur_d < 0 | fish$SpawnSitediff < -2){fish$SpawnSiteSuccess <- 0} #within 2 km spatial resolution

        #SIMULATION OF ENERGY USE DURING HOLDING
        #fish holds on the spawning grounds at maintenance metabolism until it spawns
        #spawn date drawn from the redd-survey distribution, truncated to on or after arrival;
        #holding energy averaged over nspawn draws
        #Pm_hold is the sum of one hourly Pm per day held (x24 converts to a daily total in script 10)
        A <- round(fish$Arrival_J)
        reddrow <- which.min(abs(data_ID$RKM - fish$RKM))
        cumPm <- cumsum(df_Pm[reddrow,]) #cumulative hourly Pm, by day, at the spawning grounds
        pool <- Spawning[Spawning >= A & Spawning <= ncol(df_t)]
        if(length(pool) == 0){
          fish$Pm_hold <- 0 #arrived after the last possible spawn date
        }else{
          Sdraw <- pool[((1:nspawn - 1) %% length(pool)) + 1] #resample the eligible spawn dates
          Sidx <- pmin(pmax(round(Sdraw), 1), ncol(df_t))   #spawn day
          fish$Pm_hold <- mean(cumPm[Sidx] - ifelse(A > 1, cumPm[A - 1], 0)) #mgO2/kg, sum of hourly Pm from arrival to spawning
        }
        fish$Year <- year
        fish$duty <- duty
        MH <- rbind(MH, fish)
      }
    }
    if(names(scenarios)[s] == "Figure"){
      write.csv(MH, "Output/Sims/MH_Figure.csv", row.names = FALSE)
    }else{
      write.csv(MH, paste0("Output/Sims/MH_", names(scenarios)[s], "_level", lev, ".csv"), row.names = FALSE)
    }
    print(c(names(scenarios)[s], lev, round(duty, 3), round(mean(MH$COST)*1.358442*(10^-5), 2)))
  }
}
#####

print(Sys.time() - starttime)

##END#####
