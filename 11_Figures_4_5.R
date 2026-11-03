#Figs for energy use and spawning success by date of freshwater entry (passing Rio Vista, RKM 0)
#Figure 4: energy required by stage, stacked, against initial energy density
#Figure 5: predicted spawning success, coloured by temperature (A) and velocity (B) experienced
#One simulated fish per entry date, at the fitted mean moving fraction (unweighted by entry timing)
#Results, Figure 4, Figure 5; must run script 9 first
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#####packages required
#####
library(ggplot2)
library(ggtext)
library(cowplot)
library(colorRamps)
library(scales)
library(dplyr)
library(truncnorm)
#####

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
if(!file.exists("Output/Sims/MH_Figure.csv")) stop("Run 9_Emm_Eh_sims.R first")

j = 1.358442*(10^-5) #MJ/mgO2 #to convert from mgO2 to MJ

#read in data#####
#one simulated fish per entry date (DOY 180-365) per year, at the fitted mean moving fraction
MH_all <- read.csv("Output/Sims/MH_Figure.csv")
dim(MH_all) #372
n <- nrow(MH_all)

#observed dates tracked fish passed Rio Vista
RV <- read.csv("Data/Telemetry_RioVista_passage_2022_2023.csv")
RV$jday <- as.numeric(format(as.POSIXct(RV$RV_passage_date_UTC, tz="UTC"), "%j"))
summarise(group_by(RV, Year), first = min(jday), last = max(jday), n = n())
#####

#Energy budget for each simulated fish#####
#same parameters as script 10 (Table 2); 10,000 replicates for each entry date
Ei_mean <- 9.5
Ei_sd <- 0.4
n_reps <- 10000

Emm <- MH_all$COST*j        #mainstem river migration
Eh  <- MH_all$Pm_hold*24*j  #holding

set.seed(1)
Ei_sum <- Embd_sum <- Eg_sum <- Ed_sum <- Es_sum <- SpnSuccess_sum <- rep(0, n)
for(i in 1:n_reps){
  Ei        <- sample(rnorm(1000, mean=Ei_mean, sd=Ei_sd),                size=n, replace = TRUE)
  pct_alloc <- sample(rnorm(1000, mean=0.0783, sd=0.0028),                size=n, replace = TRUE)
  Embd      <- sample(rlnorm(1000, meanlog=0.028, sdlog=0.590),           size=n, replace = TRUE)
  Eg        <- pct_alloc*Ei
  Ed        <- sample(rnorm(1000, mean=4.4, sd=0.7),                      size=n, replace = TRUE)
  Es        <- sample(rtruncnorm(1000, a=0, mean=1.30, sd=0.84),          size=n, replace = TRUE)
  EnUse <- Ei - (Embd+Emm+Eg+Eh+Es+Ed)
  SpnSuccess_val <- ifelse(EnUse >= 0 & MH_all$SpawnSiteSuccess == 1, 1, 0)
  Embd_sum <- Embd_sum + Embd
  Eg_sum   <- Eg_sum + Eg
  Ed_sum   <- Ed_sum + Ed
  Es_sum   <- Es_sum + Es
  SpnSuccess_sum <- SpnSuccess_sum + SpnSuccess_val
}
#mean of each term for each simulated fish
MH_all$E_D   <- Ed_sum/n_reps
MH_all$E_Mbd <- Embd_sum/n_reps
MH_all$E_Mr  <- Emm
MH_all$E_H   <- Eh
MH_all$E_G   <- Eg_sum/n_reps
MH_all$E_S   <- Es_sum/n_reps
MH_all$SpnSuccess_pct <- SpnSuccess_sum/n_reps*100

#date for plotting; all years shown on a common (non-leap) calendar
MH_all$pdate <- as.Date(MH_all$EntryDate_J - 1, origin = "2025-01-01")
XLIM <- as.Date(c("2025-07-28", "2025-12-03"))
#####

#FIGURE 4#####
###cumulative energy budget, fish that reached the spawning grounds; E_D at bottom, then E_Mbd, E_Mr, E_H, E_G, E_S
stage <- subset(MH_all, SpawnSiteSuccess == 1 & pdate >= XLIM[1] & pdate <= XLIM[2])
comps <- c("E_D", "E_Mbd", "E_Mr", "E_H", "E_G", "E_S")
cum <- t(apply(as.matrix(stage[, comps]), 1, cumsum))
long <- data.frame()
for(k in 1:length(comps)){
  if(k == 1){ymin <- 0}else{ymin <- cum[, k - 1]}
  long <- rbind(long, data.frame(Year = stage$Year, pdate = stage$pdate, ymin = ymin, ymax = cum[, k],
                                 component = comps[k]))
}
long$component <- factor(long$component, levels = comps)

#where the stack peaks, MJ/kg
summarise(group_by(data.frame(Year = stage$Year, top = cum[, 6]), Year), max_top = max(top))

###entry dates on which the simulated fish could not advance at some point (held at maintenance)
halt <- subset(MH_all, Days_halted > 0)
summarise(group_by(halt, Year), last_halted_entry = max(EntryDate_J), max_days_halted = max(Days_halted))
#no entry date within the observed range of Rio Vista passage was affected (both should be TRUE)
max(halt$EntryDate_J[halt$Year == 2022]) < min(RV$jday[RV$Year == 2022])
max(halt$EntryDate_J[halt$Year == 2023]) < min(RV$jday[RV$Year == 2023])
#one tile per halted entry date within the plotted range
halt <- subset(halt, pdate > XLIM[1] & pdate < XLIM[2])

###observed Rio Vista passage, shaded bar along the base
od <- summarise(group_by(RV, Year), xmin = as.Date(min(jday) - 1, origin = "2025-01-01"),
                xmax = as.Date(max(jday) - 1, origin = "2025-01-01"))

###legend labels and colours #colour blind friendly palette
bar_lab  <- "Observed RV detections"
halt_lab <- "Model: fish unable to advance"
labs <- c(E_D   = "**E<sub>D</sub>** &nbsp;energy floor (death)",
          E_Mbd = "**E<sub>Mbd</sub>** &nbsp;bay/delta migration",
          E_Mr  = "**E<sub>Mr</sub>** &nbsp;river migration",
          E_H   = "**E<sub>H</sub>** &nbsp;holding",
          E_G   = "**E<sub>G</sub>** &nbsp;gonad",
          E_S   = "**E<sub>S</sub>** &nbsp;spawning")
pal <- c(E_D = "#0072B2", E_Mbd = "#56B4E9", E_Mr = "#009E73", E_H = "#F0E442", E_G = "#E69F00", E_S = "#CC79A7")
fill_lv  <- c(bar_lab, halt_lab, comps)
fill_val <- c("grey35", "#D55E00", pal)
names(fill_val) <- fill_lv
fill_lab <- c(bar_lab, halt_lab, labs)
names(fill_lab) <- fill_lv
od$component   <- bar_lab
halt$component <- halt_lab

#initial energy density mean and +/- 2 SD, one set of lines per panel
ei_df <- data.frame(Year = rep(c(2022, 2023), each = 3),
                    yint = rep(c(Ei_mean - 2*Ei_sd, Ei_mean, Ei_mean + 2*Ei_sd), 2),
                    key  = rep(c("±2 SD", "mean", "±2 SD"), 2))

Stageplot <-
  ggplot() +
  #E_I +/- 2 SD band, under the ribbons
  annotate("rect", xmin = XLIM[1], xmax = XLIM[2], ymin = Ei_mean - 2*Ei_sd, ymax = Ei_mean + 2*Ei_sd,
           fill = "grey55", alpha = 0.18) +
  geom_ribbon(data = long, aes(x = pdate, ymin = ymin, ymax = ymax, fill = component),
              colour = "white", linewidth = 0.12) +
  geom_hline(data = ei_df, aes(yintercept = yint, linetype = key), colour = "grey15", linewidth = 0.35) +
  #observed detections and halted entry dates, in the margin below zero
  geom_rect(data = od, aes(xmin = xmin, xmax = xmax, ymin = -0.80, ymax = -0.30, fill = component), colour = NA) +
  geom_rect(data = halt, aes(xmin = pdate - 0.5, xmax = pdate + 0.5, ymin = -1.45, ymax = -0.95, fill = component),
            colour = NA) +
  scale_fill_manual(values = fill_val, labels = fill_lab, limits = fill_lv, breaks = fill_lv, drop = FALSE,
                    name = "Energy budget term",
                    guide = guide_legend(reverse = TRUE, order = 1, byrow = TRUE)) +
  scale_linetype_manual(values = c(mean = "solid", `±2 SD` = "22"), limits = c("mean", "±2 SD"),
                        name = "Initial energy density (E<sub>I</sub>)",
                        guide = guide_legend(order = 2, override.aes = list(colour = "grey15"))) +
  scale_x_date(limits = XLIM, breaks = seq(as.Date("2025-08-01"), as.Date("2025-12-01"), by = "1 month"),
               minor_breaks = NULL, labels = date_format("%b"), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 21, 3), expand = expansion(mult = c(0.01, 0.06))) +
  labs(x = "Date of passing Rio Vista (RKM 0)", y = "Energy required (MJ kg<sup>−1</sup>)") +
  #clip with coord_cartesian, not scale limits, which would remove rows from the stacked ribbons
  coord_cartesian(ylim = c(-1.55, 19)) +
  facet_wrap(~Year, ncol = 2) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border       = element_rect(colour = "grey25", linewidth = 0.4),
        panel.spacing      = unit(1.1, "lines"),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 12, margin = margin(b = 5)),
        axis.ticks         = element_line(colour = "grey25", linewidth = 0.3),
        axis.ticks.length  = unit(3, "pt"),
        axis.title.y       = element_markdown(margin = margin(r = 8)),
        axis.title.x       = element_text(margin = margin(t = 6)),
        axis.text          = element_text(colour = "grey20"),
        legend.text        = element_markdown(size = 9.5, lineheight = 1.1),
        legend.title       = element_markdown(size = 9.5, face = "bold"),
        legend.key.size    = unit(0.95, "lines"),
        legend.key         = element_blank(),
        legend.spacing.y   = unit(4, "pt"),
        legend.margin      = margin(l = 2),
        plot.margin        = margin(6, 8, 4, 6))

ggsave("Output/Figure4_EnergyUsebyStage.jpeg", Stageplot, width=10, height=5.2, units='in', dpi = 400)
#####

#FIGURE 5#####
###spawning success by entry date, coloured by mean temperature and velocity over the whole
###migration and holding period
succ <- subset(MH_all, pdate >= XLIM[1] & pdate <= XLIM[2])
succ$Year <- factor(succ$Year)
succ$Temp <- succ$Temp_mean
succ$Velocity <- succ$Uc_mean

#peak predicted spawning success, % and entry date
max(succ$SpnSuccess_pct)
succ$pdate[which.max(succ$SpnSuccess_pct)]

theme_set(theme_cowplot())
fig_temp <-
  ggplot(data = succ) +
  #thin black outline so overlapping points separate; shape distinguishes years
  geom_point(aes(x = pdate, y = SpnSuccess_pct, fill = Temp, shape = Year),
             size = 3, colour = "black", stroke = 0.30) +
  scale_fill_gradientn(colours = matlab.like2(7), limits = c(5, 20), na.value = "darkred", name = "Temp (C)") +
  scale_shape_manual(values = c(21, 24)) +
  guides(shape = guide_legend(override.aes = list(fill = "grey70"))) +
  scale_x_date(limits = XLIM, date_breaks = "1 month", labels = date_format("%b")) +
  xlab("Date of passing Rio Vista (RKM 0)") + ylab("Predicted spawning success (%)")

fig_vel <-
  ggplot(data = succ) +
  geom_point(aes(x = pdate, y = SpnSuccess_pct, fill = Velocity, shape = Year),
             size = 3, colour = "black", stroke = 0.30) +
  scale_fill_gradientn(colours = matlab.like2(7), limits = c(0, 0.4), na.value = "darkred", name = "Uc (m/s)") +
  scale_shape_manual(values = c(21, 24)) +
  guides(shape = guide_legend(override.aes = list(fill = "grey70"))) +
  scale_x_date(limits = XLIM, date_breaks = "1 month", labels = date_format("%b")) +
  xlab("Date of passing Rio Vista (RKM 0)") + ylab("Predicted spawning success (%)")

fig_succ <- plot_grid(fig_temp, fig_vel, ncol=2, labels = "AUTO")
ggsave("Output/Figure5_Entryday_SpnSucc_byTempUc.jpeg", fig_succ, width=10, height=4.5, units='in', dpi = 300)
#####

###END###
