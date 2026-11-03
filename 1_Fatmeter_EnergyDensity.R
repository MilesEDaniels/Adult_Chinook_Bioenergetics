#Energy density of adult fall-run Chinook salmon by migration stage
#Converts fatmeter scores to somatic energy density; summary statistics by stage and year;
#ANOVA and Tukey HSD tests between stages and between years within each stage
#Results, Table 1, Figure 3; EI and ED in Table 2
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#####packages required
#####
library(ggplot2)
library(dplyr)
library(ggtext)
library(multcompView)
#####

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
dir.create("Output", showWarnings = FALSE)

#read data, convert fatmeter to somatic energy density#####
#fall-run fish sampled 2021-2024; one row per fish
data <- read.csv("Data/Fatmeter_EnergyDensity_FallRun.csv")
dim(data) #3954, number of fish samples

#average of the two fatmeter scans
data$Average_1_2 <- rowMeans(data[, c("S1","S2")])

##FATMETER:ENERGY function
#Crossin & Hinch (2005); see Supplement 2 for comparison with other functions
crossin.hinch.function <- function(x) 1.7338*log(x) + 5.4886

#apply lipid:energy equation to our fatmeter data #MJ/kg
data$crossin.hinch <- crossin.hinch.function(data$Average_1_2)
#####

#Define stages#####
#MigrationReady (stage 1): fish captured in the ocean, presumed to be entering freshwater soon
#Migration (stage 2): fish migrating, but not yet on the spawning grounds
#SpawnReady (stage 3): at hatchery, fish that were spawn-ready but sacrificed prior to spawning
#Spawning (stage 4): at hatchery, fish that were killed immediately after spawning
#Postspawned (stage 5): post-spawned carcasses that completed the full spawning cycle
#Prespawnmort (stage 6): fish found dead prior to spawning, by "natural" causes
data$Energy_stage <- factor(data$Energy_stage,
                            levels=c("Stage_1_MigrationReady",
                                     "Stage_2_Migration",
                                     "Stage_3_SpawnReady",
                                     "Stage_4_Spawning",
                                     "Stage_5_Postspawned",
                                     "Stage_6_Prespawnmort"))

###Stanislaus River stage 2 fish
##these fish have suspiciously low fatmeter scores for migrating fish,
##including some below reported energy at death; based on consensus, redefine as stage 3
##Table 1 and all statistical tests include these fish as stage 3;
##the published Figure 3 plots stage 3 without them (see FIGURE 3 below)
data$STAN_reclass <- data$Energy_stage == "Stage_2_Migration" & data$Location_abrev == "STAN"
sum(data$STAN_reclass) #575
data[data$STAN_reclass == TRUE, "Energy_stage"] <- "Stage_3_SpawnReady"
summary(data$Energy_stage)
#####

#TABLE 1#####
#number of fish, fork length and energy density by stage, all years
fallrun_table <- summarise(group_by(data, Energy_stage),
                           n = n(),
                           FL = mean(FL_mm, na.rm=TRUE), FL_SD = sd(FL_mm, na.rm=TRUE),
                           Fatmeter = mean(Average_1_2), Fatmeter_SD = sd(Average_1_2),
                           EnergyDensity = mean(crossin.hinch), EnergyDensity_SD = sd(crossin.hinch), .groups = "drop")
as.data.frame(fallrun_table)
sum(fallrun_table$n) #3954
write.csv(fallrun_table, "Output/Table1_EnergyDensity_byStage.csv", row.names = FALSE)

#also by stage and year
fallrun_table_yr <- summarise(group_by(data, Energy_stage, Study_Year),
                              n = n(),
                              FL = mean(FL_mm, na.rm=TRUE), FL_SD = sd(FL_mm, na.rm=TRUE),
                              EnergyDensity = mean(crossin.hinch), EnergyDensity_SD = sd(crossin.hinch), .groups = "drop")
write.csv(fallrun_table_yr, "Output/EnergyDensity_byStage_Year.csv", row.names = FALSE)

###TABLE 2 parameters from this script
#EI, energy at the initiation of migration = stage 1
(Ei_mean <- fallrun_table$EnergyDensity[fallrun_table$Energy_stage == "Stage_1_MigrationReady"]) #9.5
(Ei_sd   <- fallrun_table$EnergyDensity_SD[fallrun_table$Energy_stage == "Stage_1_MigrationReady"]) #0.4
#ED, energy at death = post-spawned carcasses (stage 5)
(Ed_mean <- fallrun_table$EnergyDensity[fallrun_table$Energy_stage == "Stage_5_Postspawned"]) #4.4
(Ed_sd   <- fallrun_table$EnergyDensity_SD[fallrun_table$Energy_stage == "Stage_5_Postspawned"]) #0.7
#####

#ANOVA and Tukey HSD#####
###between stages, years pooled
aov_res <- aov(crossin.hinch ~ Energy_stage, data = data)
summary(aov_res)
tukey_stage <- as.data.frame(TukeyHSD(aov_res)$Energy_stage)
tukey_stage$comparison <- rownames(tukey_stage)
tukey_stage[, c("diff", "p adj")]
write.csv(tukey_stage, "Output/TukeyHSD_stage_results.csv", row.names = FALSE)

###between stage-year groups
data$stageYr <- paste0(data$Energy_stage, "_", data$Study_Year)
aov_res <- aov(crossin.hinch ~ stageYr, data = data)
summary(aov_res)
tukey_stageYr <- as.data.frame(TukeyHSD(aov_res)$stageYr)
tukey_stageYr$comparison <- rownames(tukey_stageYr)
#pull the two groups out of each comparison, e.g. "Stage_3_SpawnReady_2023-Stage_3_SpawnReady_2022"
tukey_stageYr$group1 <- sub("-.*", "", tukey_stageYr$comparison)
tukey_stageYr$group2 <- sub(".*-", "", tukey_stageYr$comparison)
tukey_stageYr$stage1 <- substr(tukey_stageYr$group1, 1, nchar(tukey_stageYr$group1) - 5)
tukey_stageYr$stage2 <- substr(tukey_stageYr$group2, 1, nchar(tukey_stageYr$group2) - 5)
tukey_stageYr$Year1  <- substr(tukey_stageYr$group1, nchar(tukey_stageYr$group1) - 3, nchar(tukey_stageYr$group1))
tukey_stageYr$Year2  <- substr(tukey_stageYr$group2, nchar(tukey_stageYr$group2) - 3, nchar(tukey_stageYr$group2))
write.csv(tukey_stageYr, "Output/TukeyHSD_stage_YR_results.csv", row.names = FALSE)
#####

#Significance letters for Figure 3#####
#stage labels for the figure
stage_labels <- c(Stage_1_MigrationReady = "Migration-ready", Stage_2_Migration = "Migration",
                  Stage_3_SpawnReady = "Spawn-ready",       Stage_4_Spawning  = "Spawning",
                  Stage_5_Postspawned = "Post-spawned",     Stage_6_Prespawnmort = "Pre-spawn\nmortality")
data$stage <- factor(stage_labels[as.character(data$Energy_stage)], levels = unname(stage_labels))
data$Year <- factor(data$Study_Year)

###between stages (bold capitals above each box)
#pre-spawn mortalities (stage 6, n = 5) not evaluated because of small sample size
#multcompLetters splits names on "-", so stages are coded S1-S5 for the letters
pv <- tukey_stage$`p adj`
names(pv) <- tukey_stage$comparison
pv <- pv[!grepl("Stage_6", names(pv))]
for(k in 1:5){
  names(pv) <- gsub(levels(data$Energy_stage)[k], paste0("S", k), names(pv))
}
between <- multcompLetters(pv)$Letters
between <- between[paste0("S", 1:5)]
#relabel so that letters read a, b, c... from left to right; groups sharing a letter still share one
between <- letters[match(between, unique(between))]
names(between) <- unname(stage_labels[1:5])
between #a, b, c, c, d

###between years within each stage (small italics beneath the points)
#shown only for stages in which at least two years differed
within <- data.frame()
for(k in 1:5){
  st <- levels(data$Energy_stage)[k]
  tk <- subset(tukey_stageYr, stage1 == st & stage2 == st)
  if(nrow(tk) == 0) next #only one year sampled
  pv <- tk$`p adj`
  names(pv) <- paste0("Y", tk$Year1, "-Y", tk$Year2)
  lett <- multcompLetters(pv)$Letters
  lett <- lett[sort(names(lett))] #year order
  lett <- letters[match(lett, unique(lett))]
  if(length(unique(lett)) > 1){
    within <- rbind(within, data.frame(stage = unname(stage_labels[k]),
                                       Year = sort(unique(c(tk$Year1, tk$Year2))),
                                       lab = lett))
  }
}
within #Migration, Spawn-ready, Spawning
#####

#FIGURE 3#####
###the published figure plots stage 3 WITHOUT the 575 reclassified Stanislaus fish (n = 1,008),
###whereas Table 1 and the significance letters above include them (n = 1,583)
figdata <- subset(data, STAN_reclass == FALSE)
table(figdata$stage)

###letter positions
#between-stage letters just above the 99th percentile of each stage
top <- summarise(group_by(figdata, stage), y = quantile(crossin.hinch, 0.99))
top$lab <- toupper(between[as.character(top$stage)])
top <- subset(top, !is.na(lab))
#within-stage letters beneath the lowest point in each stage, one baseline per stage
#x position follows the dodge of the points by year
DODGE <- 0.78
bot <- merge(within, summarise(group_by(figdata, stage), y = min(crossin.hinch)), by = "stage")
bot$x <- match(bot$stage, levels(data$stage)) +
         (match(bot$Year, levels(data$Year)) - (nlevels(data$Year) + 1)/2) * (DODGE/nlevels(data$Year))

#colour blind friendly palette (Okabe-Ito)
year_cols <- c("2021"="#E69F00", "2022"="#009E73", "2023"="#0072B2", "2024"="#CC79A7")

set.seed(700)
fig <-
  ggplot(figdata, aes(stage, crossin.hinch)) +
  #points jittered slightly in both axes to reduce overplotting; fatmeter-derived energy densities
  #cluster on repeated values. 0.035 MJ/kg is ~0.5% of the plotted range
  geom_point(aes(colour = Year), alpha = 0.55, shape = 16, size = 1.15,
             position = position_jitterdodge(jitter.width = 0.42, jitter.height = 0.035,
                                             dodge.width = DODGE)) +
  geom_boxplot(fill = NA, colour = "grey15", width = 0.52, linewidth = 0.45, outlier.shape = NA) +
  geom_text(data = top, aes(y = y, label = lab), vjust = -0.55,
            fontface = "bold", size = 3.9, colour = "grey10") +
  geom_text(data = bot, aes(x = x, y = y, label = lab), vjust = 1.9,
            fontface = "italic", size = 2.7, colour = "grey35") +
  scale_colour_manual(values = year_cols, name = "Year") +
  scale_y_continuous(breaks = seq(3, 10, 1), expand = expansion(mult = c(0.08, 0.06))) +
  labs(x = "Life stage", y = "Energy density (MJ kg<sup>−1</sup>)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.4, shape = 16))) +
  theme_bw(base_size = 11) +
  theme(axis.title.y = element_markdown(size = 10.5),
        axis.title.x = element_text(size = 10.5, margin = margin(t = 6)),
        axis.text.x  = element_text(size = 9.5, colour = "grey15", lineheight = 0.95),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey30", linewidth = 0.4),
        legend.position = "right", legend.key.size = unit(11, "pt"),
        legend.title = element_text(size = 9.5, face = "bold"),
        legend.text = element_text(size = 9),
        plot.margin = margin(6, 6, 4, 6))

ggsave("Output/Figure3_EnergyDensity_byStage.jpeg", fig, width=7.6, height=4.4, units='in', dpi = 400)
#####

###END###
