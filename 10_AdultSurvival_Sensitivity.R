#Spawning success of adults that pass Rio Vista, and one-at-a-time sensitivity analysis
#Energy budget (Eq. 1): a fish spawns successfully if it reaches the spawning grounds and
#EI - (EMbd + EMr + EG + EH + ES + ED) >= 0
#EMr and EH come from the simulations (script 9); the other terms are sampled from their distributions (Table 2)
#Simulated fish are weighted by observed entry timing, and survival is averaged over the nine levels
#of the moving fraction
#Results, Table 4; must run script 9 first
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#####packages required
#####
library(truncnorm)
#####

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
if(!file.exists("Output/Sims/MH_Baseline_level1.csv")) stop("Run 9_Emm_Eh_sims.R first")

starttime <- Sys.time()

j = 1.358442*(10^-5) #MJ/mgO2 #to convert from mgO2 to MJ
n_reps <- 1000 #Monte Carlo replicates

#Energy budget parameters, baseline (Table 2)#####
#all values are somatic energy density, MJ/kg
#EI: initial energy, normal, stage 1 fish (script 1)
#EG: gonad allocation, normal proportion of EI (script 3)
#EMbd: bay/delta migration, lognormal, tidal surfing (script 5)
#ED: energy at death, normal, post-spawned carcasses (script 1)
#ES: spawning cost, normal truncated at 0; mean 1.30, sd 0.84 (script 4, Table 2)
baseline <- data.frame(Ei_mean = 9.5, Ei_sd = 0.4, pct_mean = 0.0783, pct_sd = 0.0028,
                       Embd_meanlog = 0.028, Embd_sdlog = 0.590, Ed_mean = 4.4, Ed_sd = 0.7,
                       Es_mean = 1.30, Es_sd = 0.84)

#Timing of freshwater entry (date of passing Rio Vista)#####
#normal distribution pooled across the 2022 and 2023 tracked fish, truncated at the earliest observed
#entry (script 5); the same for both years, so that differences between years reflect river conditions
Entry_mean <- 259.5 #Sep 16
Entry_sd   <- 16.6
Entry_min  <- 234   #Aug 22
n_entry    <- 10000 #simulated fish per year
#####

#Model runs: baseline, then one factor varied at a time#####
#sims: which set of simulations (script 9); the remaining columns replace one baseline parameter
runs <- data.frame(
  run    = c("0", "1", "2a", "2b", "3a", "3b", "4a", "4b", "5a", "5b", "6a", "6b", "7a", "7b", "8"),
  factor = c("(baseline)", "Encounter velocity", "Pejus temperature", "Pejus temperature",
             "Moving fraction", "Moving fraction", "Energy floor E_D", "Energy floor E_D",
             "Gonad allocation", "Gonad allocation", "Initial energy E_I", "Initial energy E_I",
             "Spawning cost E_S", "Spawning cost E_S", "Bay/delta E_Mbd"),
  level  = c("10th percentile velocity, fitted moving fraction", "mean velocity",
             "-1 SD (15.74 C)", "+1 SD (21.66 C)", "mean -25%", "mean +25%",
             "3.5 MJ/kg", "4.8 MJ/kg", "female 20%", "male 4%", "9.0 MJ/kg", "10.0 MJ/kg",
             "1.0 MJ/kg", "2.0 MJ/kg", "active swimming, mean 2.2 MJ/kg"),
  sims   = c("Baseline", "VEL_MEAN", "PejusMinus1SD", "PejusPlus1SD", "MovingFraction_m25", "MovingFraction_p25",
             rep("Baseline", 9)),
  Ed_mean      = c(rep(4.4, 6), 3.5, 4.8, rep(4.4, 7)),
  pct_mean     = c(rep(0.0783, 8), 0.20, 0.04, rep(0.0783, 5)),
  Ei_mean      = c(rep(9.5, 10), 9.0, 10.0, rep(9.5, 3)),
  Es_mean      = c(rep(1.30, 12), 1.0, 2.0, 1.30),
  Embd_meanlog = c(rep(0.028, 14), 0.614)) #lognormal with the same sdlog, centred on 2.2 MJ/kg (script 5)
runs
#####

#Spawning success#####
Table4 <- data.frame()
Survival_bylevel <- data.frame()

for(r in 1:nrow(runs)){
  print(c(runs$run[r], runs$factor[r], runs$level[r]))
  for(lev in 1:9){
    MH <- read.csv(paste0("Output/Sims/MH_", runs$sims[r], "_level", lev, ".csv"))

    ###weight simulated fish by observed entry timing
    #each simulated fish represents one entry date; resample entry dates from the entry distribution
    set.seed(975)
    MH_all <- data.frame()
    for(year in c(2022, 2023)){
      MH_yr <- subset(MH, Year == year)
      entry <- as.integer(rtruncnorm(n_entry, a = Entry_min, mean = Entry_mean, sd = Entry_sd))
      MH_all <- rbind(MH_all, MH_yr[match(entry, MH_yr$EntryDate_J),])
    }
    MH_all <- subset(MH_all, !is.na(EntryDate_J)) #entry dates outside the simulated range, if any

    #EMr and EH from adult simulation models
    Emm <- MH_all$COST*j        #mainstem river migration
    Eh  <- MH_all$Pm_hold*24*j  #holding = sum of daily maintenance over days held
    n <- nrow(MH_all)

    ###Monte Carlo: sample other parameters, and count fish with energy left over (EnUse >= 0)
    #that also reached the spawning grounds
    set.seed(1)
    SpnSuccess_sum <- rep(0, n)
    for(i in 1:n_reps){
      Ei        <- sample(rnorm(1000, mean=runs$Ei_mean[r], sd=baseline$Ei_sd),       size=n, replace = TRUE)
      pct_alloc <- sample(rnorm(1000, mean=runs$pct_mean[r], sd=baseline$pct_sd),     size=n, replace = TRUE)
      Embd      <- sample(rlnorm(1000, meanlog=runs$Embd_meanlog[r], sdlog=baseline$Embd_sdlog), size=n, replace = TRUE)
      Eg        <- pct_alloc*Ei #gonad formation is a proportion of initial somatic energy
      Ed        <- sample(rnorm(1000, mean=runs$Ed_mean[r], sd=baseline$Ed_sd),       size=n, replace = TRUE)
      Es        <- sample(rtruncnorm(1000, a=0, mean=runs$Es_mean[r], sd=baseline$Es_sd), size=n, replace = TRUE)

      EnUse <- Ei - (Embd+Emm+Eg+Eh+Es+Ed) #energy remaining
      SpnSuccess_val <- ifelse(EnUse >= 0 & MH_all$SpawnSiteSuccess == 1, 1, 0)
      SpnSuccess_sum <- SpnSuccess_sum + SpnSuccess_val
    }
    SpnSuccess <- SpnSuccess_sum/n_reps #probability of spawning success for each simulated fish

    Survival_bylevel <- rbind(Survival_bylevel,
                              data.frame(run = runs$run[r], level = lev, duty = MH_all$duty[1],
                                         overall = mean(SpnSuccess)*100,
                                         y2022 = mean(SpnSuccess[MH_all$Year == 2022])*100,
                                         y2023 = mean(SpnSuccess[MH_all$Year == 2023])*100,
                                         Emr_mean = mean(Emm), Eh_mean = mean(Eh)))
  }
  #spawning success = mean across the nine levels of the moving fraction
  surv <- subset(Survival_bylevel, run == runs$run[r])
  Table4 <- rbind(Table4, data.frame(run = runs$run[r], factor = runs$factor[r], level = runs$level[r],
                                     overall = mean(surv$overall), y2022 = mean(surv$y2022), y2023 = mean(surv$y2023)))
  print(Table4[nrow(Table4),])
}
Table4$delta <- Table4$overall - Table4$overall[1]
#####

#RESULTS#####
###baseline spawning success, % of adults that passed Rio Vista
Table4[1,] #overall, 2022, 2023
Table4$y2023[1]/Table4$y2022[1] #year contrast, 2023:2022

###mean energy used during river migration and holding, entry-weighted, MJ/kg
colMeans(subset(Survival_bylevel, run == "0")[c("Emr_mean", "Eh_mean")])

###TABLE 4: one-at-a-time sensitivity, range of predicted spawning success for each factor
Table4
write.csv(Table4, "Output/Table4_Sensitivity_OAT.csv", row.names = FALSE)
write.csv(Survival_bylevel, "Output/Survival_bylevel.csv", row.names = FALSE)
#####

print(Sys.time() - starttime)

###END###
