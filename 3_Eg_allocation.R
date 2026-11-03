#Percent of initial somatic energy that females allocate to gonad formation (EG)
#Klamath River fall-run females (Hearsey 2011) used as a proxy for Central Valley fall-run females
#Results, Table 2, Supplement 4
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
dir.create("Output", showWarnings = FALSE)

#read data for model from Klamath/Trinity R#####
data <- read.csv("Data/Klamath_Hearsey2011_FL_Mass_Gonad.csv")
dim(data) #704

data$Gonad_mass_kg <- data$Gonad_mass_g / 1000
data$Gonad_Energy_MJ <- data$Gonad_mass_kg * 7.79 #gonad energy density 7.79 MJ/kg, Hendry et al. 1999
data$Soma_kg <- data$Mass_kg - data$Gonad_mass_kg

#subset into stages, female fall-run only
Ent <- subset(data, Stage == "Entry"       & Sex == "f" & Run == "fall")
SpR <- subset(data, Stage == "Spawn-ready" & Sex == "f" & Run == "fall")
nrow(Ent) #134
nrow(SpR) #52
#####

#linear models of ln (log e) transformed data#####
#this method used by Hendry & Berg 1999, Crossin & Hinch 2005, Bowerman et al. 2017, etc.
###gonad energy v FL
LM_Ent_GE <- lm(log(Gonad_Energy_MJ) ~ log(FL_mm), data=Ent)
summary(LM_Ent_GE) #log(Gonad_Energy_MJ) = 3.6296*log(FL_mm) - 23.2143
LM_SpR_GE <- lm(log(Gonad_Energy_MJ) ~ log(FL_mm), data=SpR)
summary(LM_SpR_GE) #log(Gonad_Energy_MJ) = 2.9546*log(FL_mm) - 17.6602

###somatic mass v FL, at freshwater entry
LM_Ent_MS <- lm(log(Soma_kg) ~ log(FL_mm), data=Ent)
summary(LM_Ent_MS)

###predict values, back-transformed, for every fork length 200-1500 mm
newdata <- data.frame(FL_mm = seq(200, 1500, 1))
Ent_SpR <- data.frame(FL_mm         = newdata$FL_mm,
                      Gon_En_Ent    = exp(predict(LM_Ent_GE, newdata = newdata)),
                      Gon_En_SpR    = exp(predict(LM_SpR_GE, newdata = newdata)),
                      Mass_soma_Ent = exp(predict(LM_Ent_MS, newdata = newdata)))

#gonad energy gained between freshwater entry and spawn-ready, MJ
Ent_SpR$Gon_En <- Ent_SpR$Gon_En_SpR - Ent_SpR$Gon_En_Ent

#somatic energy (MJ) at entry is somatic mass * somatic energy density at entry (EI)
#EI = 9.5 MJ/kg, stage 1 mean (script 1, Table 1)
Ei_mean <- 9.5
Ei_sd <- 0.4
Ent_SpR$Soma_Ent_En <- Ent_SpR$Mass_soma_Ent * Ei_mean

#percent of initial somatic energy allocated to gonad formation
Ent_SpR$Gon_alloc_pct <- Ent_SpR$Gon_En / Ent_SpR$Soma_Ent_En * 100
summary(Ent_SpR$Gon_alloc_pct)
#####

#predict % allocation for our fall-run females from their fork lengths#####
ourdata <- read.csv("Data/Fatmeter_EnergyDensity_FallRun.csv")
ourdata <- subset(ourdata, !is.na(FL_mm) & Sex_MF == "Female")
nrow(ourdata) #1756

ourdata$Gon_alloc_pct <- Ent_SpR$Gon_alloc_pct[match(as.integer(ourdata$FL_mm), Ent_SpR$FL_mm)]
summary(ourdata$Gon_alloc_pct)

###RESULTS, Table 2
(Eg_pct_mean <- mean(ourdata$Gon_alloc_pct)) #7.83 %
(Eg_pct_sd   <- sd(ourdata$Gon_alloc_pct))   #0.28 %

#EG in MJ/kg, for a fish with EI 9.5 +/- 0.4 MJ/kg
(Eg_mean <- Eg_pct_mean/100 * Ei_mean) #0.7
(Eg_sd   <- sqrt((Eg_pct_mean/100 * Ei_sd)^2 + (Eg_pct_sd/100 * Ei_mean)^2)) #0.04

Eg_out <- data.frame(Eg_pct_mean = Eg_pct_mean, Eg_pct_sd = Eg_pct_sd,
                     Eg_MJ_kg_mean = Eg_mean, Eg_MJ_kg_sd = Eg_sd, n_females = nrow(ourdata))
write.csv(Eg_out, "Output/Eg_allocation.csv", row.names = FALSE)
#####

###END###
