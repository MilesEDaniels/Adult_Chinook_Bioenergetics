#Energy density used for spawning (ES)
#Difference in somatic energy density between spawn-ready adults (stage 3) and
#post-spawned carcasses (stage 5); must run script 1 first
#Results, Table 2, Supplement 3
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
if(!file.exists("Output/Table1_EnergyDensity_byStage.csv")) stop("Run 1_Fatmeter_EnergyDensity.R first")

#read stage summary from script 1#####
fallrun_table <- read.csv("Output/Table1_EnergyDensity_byStage.csv")
fallrun_table[c("Energy_stage", "n", "EnergyDensity", "EnergyDensity_SD")]

SpawnReady_mean  <- fallrun_table$EnergyDensity[fallrun_table$Energy_stage == "Stage_3_SpawnReady"]
SpawnReady_sd    <- fallrun_table$EnergyDensity_SD[fallrun_table$Energy_stage == "Stage_3_SpawnReady"]
Postspawned_mean <- fallrun_table$EnergyDensity[fallrun_table$Energy_stage == "Stage_5_Postspawned"]
Postspawned_sd   <- fallrun_table$EnergyDensity_SD[fallrun_table$Energy_stage == "Stage_5_Postspawned"]
#####

#Estimate Es#####
#sample spawn-ready and post-spawned energy densities, take the difference,
#and keep positive values (a fish cannot gain energy by spawning)
set.seed(234)
Es_data <- data.frame(SpawnReady  = rnorm(1000000, mean = SpawnReady_mean,  sd = SpawnReady_sd),
                      Postspawned = rnorm(1000000, mean = Postspawned_mean, sd = Postspawned_sd))
Es_data$Es <- Es_data$SpawnReady - Es_data$Postspawned
Es_data <- subset(Es_data, Es > 0)

###RESULTS, Table 2
(Es_mean <- mean(Es_data$Es)) #1.30 MJ/kg
(Es_sd   <- sd(Es_data$Es))   #0.84 MJ/kg

write.csv(data.frame(Es_mean = Es_mean, Es_sd = Es_sd), "Output/Es_parameter.csv", row.names = FALSE)
#####

###END###
