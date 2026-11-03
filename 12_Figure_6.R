#Fig for the one-at-a-time sensitivity analysis of predicted spawning success (tornado plot)
#Each factor varied alone between the levels shown, all others held at baseline (Table 4)
#Results, Figure 6; must run script 10 first
#Miles E. Daniels, miles.daniels@ucsc.edu; Alyssa FitzGerald, alyssa.fitzgerald@noaa.gov
#08-Oct-2026
#Adapted from FitzGerald, A.M., & Martin, B.M. (2022). Conservation Physiology

rm(list=ls())
gc()

#####packages required
#####
library(ggplot2)
library(ggtext)
#####

#set working directory to the Final_Code_Data folder, e.g.
#setwd("C:/Users/yourname/Documents/Final_Code_Data")
if(!dir.exists("Data")) stop("Set the working directory to the Final_Code_Data folder")
if(!file.exists("Output/Table4_Sensitivity_OAT.csv")) stop("Run 10_AdultSurvival_Sensitivity.R first")

#read in data#####
Table4 <- read.csv("Output/Table4_Sensitivity_OAT.csv", colClasses = c(run = "character"))
Table4
base <- Table4$overall[Table4$run == "0"] #baseline spawning success, %
#####

#factors and the two levels compared for each#####
#encounter velocity and bay/delta migration cost are one-sided, comparing the alternative with the baseline (run 0)
spec <- data.frame(
  factor = c("Energy floor  E<sub>D</sub>", "Gonad allocation  E<sub>G</sub>", "Pejus temperature",
             "Encounter velocity", "Initial energy  E<sub>I</sub>", "Spawning cost  E<sub>S</sub>",
             "Bay/delta  E<sub>Mbd</sub>", "Moving fraction"),
  runA   = c("4b", "5a", "2a", "1", "6a", "7b", "8", "3a"),
  labA   = c("4.8 MJ kg<sup>−1</sup>", "female 20%", "−1 SD", "mean", "9.0 MJ kg<sup>−1</sup>", "2.0 MJ kg<sup>−1</sup>", "active cruise", "−25%"),
  runB   = c("4a", "5b", "2b", "0", "6b", "7a", "0", "3b"),
  labB   = c("3.5 MJ kg<sup>−1</sup>", "male 4%", "+1 SD", "10th pct", "10.0 MJ kg<sup>−1</sup>", "1.0 MJ kg<sup>−1</sup>", "tidal surfing", "+25%"))

fac <- data.frame()
for(k in 1:nrow(spec)){
  vA <- Table4$overall[Table4$run == spec$runA[k]]
  vB <- Table4$overall[Table4$run == spec$runB[k]]
  if(vA <= vB){
    fac <- rbind(fac, data.frame(factor = spec$factor[k], lo = vA, hi = vB, lo_lab = spec$labA[k], hi_lab = spec$labB[k]))
  }else{
    fac <- rbind(fac, data.frame(factor = spec$factor[k], lo = vB, hi = vA, lo_lab = spec$labB[k], hi_lab = spec$labA[k]))
  }
}
fac$span <- fac$hi - fac$lo
fac <- fac[order(fac$span),] #ascending, so the widest bar is on top
fac$y <- 1:nrow(fac)
#the moving fraction is the only factor fitted to our own data
fac$grp <- ifelse(fac$factor == "Moving fraction", "fitted", "measured")
fac[nrow(fac):1, c("factor", "lo", "hi", "span")]
fac$grp[1] == "fitted" #TRUE if the moving fraction is the least influential factor
#####

#bars, split at the baseline into the part that lowers and the part that raises spawning success#####
h <- 0.32 #half bar height
rects <- data.frame()
for(k in 1:nrow(fac)){
  if(fac$lo[k] < base){
    rects <- rbind(rects, data.frame(y = fac$y[k], xmin = fac$lo[k], xmax = min(base, fac$hi[k]),
                                     key = paste0(fac$grp[k], "_lower")))
  }
  if(fac$hi[k] > base){
    rects <- rbind(rects, data.frame(y = fac$y[k], xmin = max(base, fac$lo[k]), xmax = fac$hi[k],
                                     key = paste0(fac$grp[k], "_raise")))
  }
}
#value labels at both ends of each bar; asterisk on the pejus -1 SD value (see footnote)
labs <- rbind(data.frame(y = fac$y, x = fac$lo, hj = 1.08,
                         lab = paste0(fac$lo_lab, "&nbsp;&nbsp;", sprintf("%.1f", fac$lo),
                                      ifelse(fac$factor == "Pejus temperature", "*", ""))),
              data.frame(y = fac$y, x = fac$hi, hj = -0.08,
                         lab = paste0(sprintf("%.1f", fac$hi), "&nbsp;&nbsp;", fac$hi_lab)))

#colour blind friendly palette (Okabe-Ito), as Figures 3-5
pal <- c(measured_lower = "#D55E00", measured_raise = "#0072B2",
         fitted_lower   = "#E8A87C", fitted_raise   = "#7FB3D9")
#####

#FIGURE 6#####
fig <-
  ggplot() +
  geom_vline(xintercept = base, linetype = "dashed", colour = "grey35", linewidth = 0.5) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = y - h, ymax = y + h, fill = key),
            colour = "white", linewidth = 0.3) +
  geom_richtext(data = labs, aes(x = x, y = y, label = lab, hjust = hj), size = 2.9, colour = "grey15",
                fill = NA, label.colour = NA, label.padding = unit(0, "pt")) +
  annotate("text", x = base, y = nrow(fac) + 0.72, label = sprintf("baseline %.1f%%", base),
           size = 3.2, colour = "grey25", fontface = "bold") +
  annotate("text", x = -8.6, y = 0.30, hjust = 0, size = 2.8, colour = "grey30",
           label = "* at 1 SD below the central pejus estimate no fish is predicted to survive; that bound is excluded by our own telemetry (see text).") +
  scale_fill_manual(values = pal, breaks = c("measured_lower", "measured_raise", "fitted_lower"),
                    labels = c("lowers spawning success", "raises spawning success",
                               "the moving fraction — the only fitted parameter, and the least influential"),
                    name = NULL) +
  scale_y_continuous(breaks = fac$y, labels = fac$factor, expand = expansion(add = c(0.95, 1.15))) +
  scale_x_continuous(limits = c(-9, 49), breaks = seq(0, 45, 5), expand = expansion(mult = 0)) +
  labs(x = "Predicted spawning success (%)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_markdown(size = 10, colour = "grey15"),
        axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey30", linewidth = 0.4),
        legend.position = "top", legend.justification = "left", legend.direction = "vertical",
        legend.key.size = unit(11, "pt"),
        legend.margin = margin(b = -4), legend.text = element_text(size = 9.5),
        plot.margin = margin(6, 10, 6, 6))

ggsave("Output/Figure6_OAT_tornado.jpeg", fig, width=9, height=5.8, units='in', dpi = 300)
#####

###END###
