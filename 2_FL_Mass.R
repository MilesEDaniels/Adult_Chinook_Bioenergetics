#Model to estimate mass from fork length
#Mass was difficult to measure for live adult Chinook salmon, so mass is estimated from fork length
#using Klamath River adults (Hearsey 2011); used to estimate bay/delta energy use (script 5)
#Methods, Supplement 1
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
dir.create("Output", showWarnings = FALSE)

#read data#####
#Klamath River fall- and spring-run adults, Hearsey (2011) Appendices F and G
data <- read.csv("Data/Klamath_Hearsey2011_FL_Mass_Gonad.csv")
dim(data) #704
#####

#linear model of ln (log e) transformed data#####
#this method used by Hendry & Berg 1999, Crossin & Hinch 2005, Bowerman et al. 2017, etc.
#stages, sexes and runs grouped
LM_FL_M <- lm(log(Mass_kg) ~ log(FL_mm), data=data)
summary(LM_FL_M) #log(Mass_kg) = 3.05055*log(FL_mm) - 18.40448

FL_Mass_coef <- data.frame(intercept = round(coef(LM_FL_M)[1], 5),
                           slope     = round(coef(LM_FL_M)[2], 5),
                           r2        = summary(LM_FL_M)$r.squared,
                           n         = nobs(LM_FL_M))
FL_Mass_coef
write.csv(FL_Mass_coef, "Output/FL_Mass_coefficients.csv", row.names = FALSE)

#example: mass of a 750 mm fish, kg
exp(FL_Mass_coef$slope*log(750) + FL_Mass_coef$intercept) #6.0
#####

###END###
