#Optimal migration swim speed based on temperature and velocity, for each formulation of critical power (Pc)
#For every temperature and water velocity, the swim speed that minimises cost of transport (COT),
#the proportion of time travelling vs recovering, and the resulting migration (ground) speed, Um
#Bioenergetics model 3 of Martin et al. 2015
#Methods, Supplement 5; used to build the energy landscapes (script 7)
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
dir.create("Output/OptSwim", showWarnings = FALSE, recursive = TRUE)

starttime <- Sys.time()

#parameters#####
M     = 5.1                 #kg                  #fish mass                                        #Martin et al. 2015
j     = 1.358442*(10^-5)    #MJ/mgO2             #to convert from mgO2 to MJ
c_0   = 65.51               #mgO2/h/kg           #Intercept for mass-specific maintenance rate     #Martin et al. 2015
b     = -0.217              #unitless            #Mass exponent for mass-specific maintenance rate #Martin et al. 2015
d     = 0.068               #C^-1                #Temperature exponent for maintenance rate        #Martin et al. 2015
c_1   = 420.6               #mgO2/h/kg/m^-gamma  #Intercept for swimming cost scalar (1 kg)        #Martin et al. 2015
beta  = -0.369              #unitless            #Mass exponent for mass-specific swim costs       #Martin et al. 2015
gamma = 2.130               #unitless            #Velocity exponent for mass-specific swim costs   #Martin et al. 2015
Nmax  = 103                 #mgO2/kg             #Anaerobic work capacity                          #Martin et al. 2015
yr    = 1.82                #unitless            #Recovery yield                                   #Martin et al. 2015

#aerobic scope of pre-spawning adult Chinook salmon, coastal fall-run population (Chilliwack River)
#Van Wert et al. 2023: temperature at maximum aerobic scope, maximum absolute aerobic scope (AAS),
#and upper pejus temperature (temperature at which AAS falls to 80% of its maximum, mean +/- 1 SD)
Topt    = 14.75             #C
AASmax  = 9.6               #mgO2/kg/min
Tpejus  = 18.70             #C
Tpej_sd = 2.96              #C
#####

#functions#####
#Pm and Ps are energy/h, mgO2/h/kg
Pm.func <- function(temp){return(c_0*(M^b)*exp(d*temp))} #Energy needed for maintenance metabolism
Ps.func <- function(Us){return(c_1*(M^beta)*Us^gamma)}   #Us is swim speed in m/s #Energy needed for swimming
#####

#Critical power (maximum aerobic capacity, Bmax or Pc) formulations#####
#"PcBase":     the temperature-dependent quadratic of FitzGerald & Martin 2022, which barely declines in
#              warm water; the base model ("base-model critical power" in Table 3)
#"PcCoastal":  aerobic scope reconstructed from Van Wert et al. 2023 as a second-order polynomial,
#              AAS(T) = AASmax - a*(T - Topt)^2, with a = 0.2*AASmax/(Tpejus - Topt)^2,
#              floored at zero; Pc = AAS + Pm. The model adopted in the paper
#"PcCoastal_pejusMinus1SD" / "PcCoastal_pejusPlus1SD": as PcCoastal with the upper pejus temperature
#              1 SD lower / higher, for the sensitivity analysis (Table 4)
Pc_names <- c("PcBase", "PcCoastal", "PcCoastal_pejusMinus1SD", "PcCoastal_pejusPlus1SD")
Pc_dpej  <- c(NA, 0, -Tpej_sd, Tpej_sd) #shift in upper pejus temperature, C
#####

#Calculate COT for a range of temperatures, velocities and swim speeds#####
temps    <- round(seq(1, 32, 0.1), 2)   #temperature, C
Uc_range <- round(seq(0.05, 3, 0.05), 2) #water velocity, m/s; velocity of 0 not used
Us_range <- round(seq(0, 5, 0.05), 2)    #swim speed, m/s; will sample out optimum Us

for(p in 1:length(Pc_names)){
  print(Pc_names[p])
  #one row per temperature and velocity
  OptSwim <- data.frame(matrix(nrow = length(temps)*length(Uc_range), ncol = 10))
  colnames(OptSwim) <- c("temp", "Uc", "Us", "Pm", "Ps", "Pc", "COT", "te", "tr", "TravProp")
  row <- 0
  for(t in 1:length(temps)){
    temp <- temps[t]
    Pm <- Pm.func(temp = temp)
    if(Pc_names[p] == "PcBase"){
      Pc <- 274.7371 + 49.0283*temp - 1.4244*temp^2
    }else{
      a <- 0.2*AASmax/((Tpejus + Pc_dpej[p]) - Topt)^2
      AAS <- max(0, AASmax - a*(temp - Topt)^2) #mgO2/kg/min
      Pc <- AAS*60 + Pm                         #mgO2/h/kg #where AAS is 0, Pc = Pm and the fish cannot migrate
    }
    for(k in 1:length(Uc_range)){
      Uc <- Uc_range[k]
      Us <- Us_range[Us_range > Uc] #remove fish swimming slower than velocity
      Ps <- Ps.func(Us = Us)
      #COT = cost of transport, energy per unit distance (mgO2/kg/km); Model 3 from Martin et al. 2015
      #COST (numerator) is energy per unit time spent swimming, including recovery (mgO2/h/kg),
      #divided by km covered per hour spent swimming [(Us-Uc)*60*60/1000]
      over <- pmax(0, Pm + Ps - Pc) #power demand above critical power, met anaerobically
      if(Pc > Pm){
        COST <- Pm + Ps + Pm*(yr*over/(Pc - Pm)) + (yr - 1)*over
      }else{
        COST <- rep(Inf, length(Us)) #no aerobic scope; no swim speed is sustainable
      }
      COT <- COST/((Us - Uc)*60*60/1000)
      best <- which.min(COT) #optimal swim speed minimises COT
      row <- row + 1
      OptSwim[row, c("temp", "Uc", "Us", "Pm", "Ps", "Pc", "COT")] <- c(temp, Uc, Us[best], Pm, Ps[best], Pc, COT[best])
    }
  }

  #predict optimal migration speeds
  OptSwim$te <- Nmax/(OptSwim$Pm + OptSwim$Ps - OptSwim$Pc) #time to exhaustion
  OptSwim$tr <- yr*Nmax/(OptSwim$Pc - OptSwim$Pm)           #time to recovery
  OptSwim[OptSwim$te <= 0, "te"] <- 0 #if te <= 0, the fish can swim at this speed indefinitely
  OptSwim$TravProp <- ifelse(OptSwim$te <= 0, 1, OptSwim$te/(OptSwim$te + OptSwim$tr)) #proportion of time travelling
  OptSwim$Um <- OptSwim$TravProp*(OptSwim$Us - OptSwim$Uc) #migration (ground) speed, m/s
  OptSwim$COT_MJ <- OptSwim$COT*j #MJ/kg/km

  print(summary(OptSwim$Um))
  #csv to read; rds keeps full numerical precision for the energy landscapes (script 7)
  write.csv(OptSwim, paste0("Output/OptSwim/OptSwim_tempvel_", M, "kg_", Pc_names[p], ".csv"), row.names = FALSE)
  saveRDS(OptSwim, paste0("Output/OptSwim/OptSwim_tempvel_", M, "kg_", Pc_names[p], ".rds"))
}
#####

print(Sys.time() - starttime)

###END###
