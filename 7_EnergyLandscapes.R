#Energy landscapes for adult migration on the mainstem Sacramento River
#Temperature, velocity, migration speed, maintenance metabolism and cost landscapes (RKM x day matrices),
#used as look-up matrices by the migration and holding simulations (scripts 8 and 9)
#River temperature and velocity from the RAFT model, Rio Vista (RKM 0) to Keswick Dam, 2 km and daily resolution
#Methods; must run script 6 first
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
if(!file.exists("Output/OptSwim/OptSwim_tempvel_5.1kg_PcCoastal.rds")) stop("Run 6_OptimalMigSwimSpeed.R first")
dir.create("Output/Landscapes", showWarnings = FALSE, recursive = TRUE)

starttime <- Sys.time()

#read in temperature and velocity data and merge into single dfs#####
#2016 is example, will cater for each year
DATE <- read.csv("Data/River/River_Time_2016.csv", header = FALSE, sep = ",") #time
DATE <- format(as.POSIXct(as.vector(DATE[,]), format="%d-%B-%Y %H:%M:%S", tz="UTC"), format="%m-%d")
DATE <- unique(DATE)
DATE[366] <- "12-31" #12-31 is not included in the time file

SPACE <- read.csv("Data/River/River_Distance_Downstream_Keswick_m_2014.csv", header = FALSE, sep = ",")
colnames(SPACE) <- "m_ds_Keswick"
#convert m to km, and reorder so that river entry (Rio Vista) is RKM 0, and highest RKM is at Keswick Dam
SPACE$m_us_mouth <- rev(SPACE$m_ds_Keswick)
SPACE$RKM <- SPACE$m_us_mouth/1000

#empty dataframes
data       <- data.frame(matrix(nrow=nrow(SPACE), ncol=0)) #temperature, C
data_vel   <- data.frame(matrix(nrow=nrow(SPACE), ncol=0)) #mean cross-sectional velocity, m/s
data_vel10 <- data.frame(matrix(nrow=nrow(SPACE), ncol=0)) #10th percentile cross-sectional velocity, m/s

#2012 is used to calibrate the moving fraction (script 8); 2022 and 2023 are the study years
#2013 is read so that 2012 fish can be followed for 50 days into the next year
for(year in c(2012, 2013, 2022, 2023)){
  temps <- read.csv(paste0("Data/River/Daily_River_Temp_C_", year, ".csv"), header = FALSE, sep = ",")
  vel   <- read.csv(paste0("Data/River/Daily_River_Velocity_m_sec_", year, ".csv"), header = FALSE, sep = ",")
  vel10 <- read.csv(paste0("Data/River/Daily_River_Velocity_m_sec_", year, "_10p.csv"), header = FALSE, sep = ",")
  #treat leap years differently
  if(year == 2012){
    DATEyr <- DATE
  }else{
    DATEyr <- DATE[!(DATE == "02-29")]
    temps <- temps[1:365] #remove last column bc NA
    vel   <- vel[1:365]
    vel10 <- vel10[1:365]
  }
  DATEyr <- paste0(DATEyr, "-", year)
  colnames(temps) <- DATEyr
  colnames(vel)   <- DATEyr
  colnames(vel10) <- DATEyr
  data       <- cbind(data, temps)
  data_vel   <- cbind(data_vel, vel)
  data_vel10 <- cbind(data_vel10, vel10)
  print(c(year, dim(temps), dim(data)))
}

###fill in NA gaps by linear interpolation
#velocities are missing at some river km (fill across RKM first), and some days (then fill across days)
data       <- as.matrix(data)
data_vel   <- as.matrix(data_vel)
data_vel10 <- as.matrix(data_vel10)
sum(is.na(data)); sum(is.na(data_vel)); sum(is.na(data_vel10))
x <- 1:nrow(data_vel)
for(k in 1:ncol(data_vel)){
  if(anyNA(data_vel[,k]) & !all(is.na(data_vel[,k]))){
    data_vel[,k] <- approx(x, data_vel[,k], xout = x, method = "linear", rule = 2)$y
  }
  if(anyNA(data_vel10[,k]) & !all(is.na(data_vel10[,k]))){
    data_vel10[,k] <- approx(x, data_vel10[,k], xout = x, method = "linear", rule = 2)$y
  }
}
x <- 1:ncol(data)
for(i in 1:nrow(data)){
  if(anyNA(data[i,])){
    data[i,] <- approx(x, data[i,], xout = x, method = "linear", rule = 2)$y
  }
  if(anyNA(data_vel[i,])){
    data_vel[i,] <- approx(x, data_vel[i,], xout = x, method = "linear", rule = 2)$y
  }
  if(anyNA(data_vel10[i,])){
    data_vel10[i,] <- approx(x, data_vel10[i,], xout = x, method = "linear", rule = 2)$y
  }
}
sum(is.na(data)); sum(is.na(data_vel)); sum(is.na(data_vel10)) #0

#reorder so that RKM 0 is at top of matrix
data       <- data[order(SPACE$RKM),]
data_vel   <- data_vel[order(SPACE$RKM),]
data_vel10 <- data_vel10[order(SPACE$RKM),]
data_ID    <- SPACE[order(SPACE$RKM),]
write.csv(data_ID, "Output/Landscapes/data_ID_Sacramento.csv", row.names = FALSE)
#####

#bioenergetics parameters and functions#####
M     = 5.1                 #kg                  #fish mass                                        #Martin et al. 2015
c_0   = 65.51               #mgO2/h/kg           #Intercept for mass-specific maintenance rate     #Martin et al. 2015
b     = -0.217              #unitless            #Mass exponent for mass-specific maintenance rate #Martin et al. 2015
d     = 0.068               #C^-1                #Temperature exponent for maintenance rate        #Martin et al. 2015
c_1   = 420.6               #mgO2/h/kg/m^-gamma  #Intercept for swimming cost scalar (1 kg)        #Martin et al. 2015
beta  = -0.369              #unitless            #Mass exponent for mass-specific swim costs       #Martin et al. 2015
gamma = 2.130               #unitless            #Velocity exponent for mass-specific swim costs   #Martin et al. 2015
yr    = 1.82                #unitless            #Recovery yield                                   #Martin et al. 2015

Pm.func <- function(temp){return(c_0*(M^b)*exp(d*temp))} #mgO2/h/kg #maintenance metabolism
Ps.func <- function(Us){return(c_1*(M^beta)*Us^gamma)}   #mgO2/h/kg #swimming, Us in m/s
#COST is energy per unit time spent swimming, mgO2/h/kg, and includes the cost of recovering from any
#anaerobic swimming; Model 3 from Martin et al. 2015
COST.func <- function(Pm, Ps, Pc){return(Pm + Ps + Pm*pmax(0, (yr*(Pm+Ps-Pc))/(Pc-Pm)) + pmax(0, (yr-1)*(Pm+Ps-Pc)))}
#####

#Landscapes for each scenario#####
#scenario = critical power formulation (script 6) and velocity metric
#10th percentile velocity: fish assumed to use the slower water at channel margins (baseline)
#mean velocity: sensitivity analysis (Table 4)
scenarios <- data.frame(name      = c("PcBase_VEL_10pct", "PcCoastal_VEL_10pct", "PcCoastal_VEL_MEAN",
                                      "PcCoastal_pejusMinus1SD_VEL_10pct", "PcCoastal_pejusPlus1SD_VEL_10pct"),
                        Pc        = c("PcBase", "PcCoastal", "PcCoastal",
                                      "PcCoastal_pejusMinus1SD", "PcCoastal_pejusPlus1SD"),
                        vel       = c("VEL_10pct", "VEL_10pct", "VEL_MEAN", "VEL_10pct", "VEL_10pct"))

#round temperature and velocity to match the OptSwim look-up tables
t_round <- round(round(data/0.10)*0.10, 2) #nearest 0.1 C
cols <- colnames(data)

for(s in 1:nrow(scenarios)){
  print(scenarios$name[s])
  OptSwim <- readRDS(paste0("Output/OptSwim/OptSwim_tempvel_", M, "kg_", scenarios$Pc[s], ".rds"))

  if(scenarios$vel[s] == "VEL_10pct"){v_data <- data_vel10}else{v_data <- data_vel}
  v_round <- round(round(v_data/0.05)*0.05, 2) #nearest 0.05 m/s
  v_round[v_round > 3] <- 3.0
  v_round[v_round == 0] <- 0.05 #velocity can't really be 0 bc then we can't calculate swim speed, so assume 0.05 m/s

  #look up the optimal swim speed, migration speed and critical power for each cell's temperature and velocity
  rown <- match(paste(t_round, v_round), paste(OptSwim$temp, OptSwim$Uc))
  if(anyNA(rown)) stop("some temperatures or velocities fall outside the OptSwim table")
  df_Um <- matrix(OptSwim$Um[rown], nrow = nrow(data)) #migration speed, m/s
  df_Us <- matrix(OptSwim$Us[rown], nrow = nrow(data)) #swim speed, m/s
  df_Pc <- matrix(OptSwim$Pc[rown], nrow = nrow(data)) #critical power, mgO2/h/kg
  df_t  <- t_round
  df_v  <- v_round
  df_Pm <- Pm.func(temp = t_round)                     #mgO2/h/kg
  df_TravProp <- matrix(OptSwim$TravProp[rown], nrow = nrow(data)) #proportion of time travelling
  #the simulation (script 9) charges energy per unit ELAPSED time, so COST is multiplied by TravProp
  #this way, energy per km equals the cost of transport (COT) wherever the fish must stop to recover
  df_COST <- COST.func(Pm = df_Pm, Ps = Ps.func(Us = df_Us), Pc = df_Pc)*df_TravProp #mgO2/h/kg
  df_COST[df_TravProp == 0] <- df_Pm[df_TravProp == 0] #no aerobic scope: the fish holds at maintenance

  #save one landscape per year: columns are day of year, plus 50 days into the next year so that
  #late entrants can still be followed (not available for 2023)
  for(year in c(2012, 2022, 2023)){
    k <- which(substr(cols, 7, 10) == as.character(year))
    if(year != 2023){k <- c(k, max(k) + 1:50)}
    landscape <- list(df_t = df_t[,k], df_v = df_v[,k], df_Um = df_Um[,k],
                      df_Pm = df_Pm[,k], df_COST = df_COST[,k], dates = cols[k])
    saveRDS(landscape, paste0("Output/Landscapes/Landscape_", scenarios$name[s], "_", year, ".rds"))
  }
  print(summary(c(df_Um)))
}
#####

print(Sys.time() - starttime)

###END###
