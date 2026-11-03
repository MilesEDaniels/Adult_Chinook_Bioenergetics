## Regenerates the model-based supplement figures from the corrected package at ES = 1.30:
## S10.1 (mean velocity budget), S13.1 (calibration), S14.1 (sensitivity by year), S15.1 (EMr split).
## Run with the working directory set to the package folder, after scripts 1-12:
##   source("Supplement_Figures/make_supplement_figures.R")
## Writes only to Supplement_Figures/; nothing in Output/ is changed.
## The original scripts are reused by patching their text, so these figures use the same code as
## the main-text figures.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
if(!file.exists("9_Emm_Eh_sims.R")) stop("Set the working directory to the package folder")
OUT <- normalizePath("Supplement_Figures")
j <- 1.358442e-5
run_patched <- function(src){ tf <- tempfile(fileext = ".R"); writeLines(src, tf)
  st <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(tf)), stdout = TRUE, stderr = TRUE)
  if(!is.null(attr(st, "status"))) {cat(tail(st, 25), sep = "\n"); stop("patched script failed")}
  invisible(st) }
fixed <- function(x, a, b){ stopifnot(any(grepl(a, x, fixed = TRUE))); gsub(a, b, x, fixed = TRUE) }

## ---- S13.1: script 8 with the simulated rates kept and its outputs redirected ----
s8 <- readLines("8_MovingFraction_Calibration.R")
s8 <- fixed(s8, "sim_rate <- sim_rate[is.finite(sim_rate)]",
            sprintf('sim_rate <- sim_rate[is.finite(sim_rate)]; saveRDS(list(sim=sim_rate, obs=tst$rate), paste0("%s/simrate_", arms[a], ".rds"))', OUT))
s8 <- fixed(s8, '"Output/Table3_MigrationRate_calibration.csv"', sprintf('"%s/Table3_check.csv"', OUT))
s8 <- fixed(s8, '"Output/MovingFraction_levels.csv"', 'tempfile()')
s8 <- fixed(s8, '"Output/MovingFraction_levels.rds"', 'tempfile()')
run_patched(s8)
chk <- read.csv(file.path(OUT, "Table3_check.csv")); pkg <- read.csv("Output/Table3_MigrationRate_calibration.csv")
stopifnot(isTRUE(all.equal(chk$KS_D, pkg$KS_D)))   # same draws as Table 3
arms <- c(PcBase = "Base-model critical power", PcCoastal = "Aerobic scope, coastal fall-run (adopted)")
sim <- bind_rows(lapply(names(arms), function(a) data.frame(arm = arms[[a]], rate = readRDS(sprintf("%s/simrate_%s.rds", OUT, a))$sim)))
obs <- bind_rows(lapply(names(arms), function(a) data.frame(arm = arms[[a]], rate = readRDS(sprintf("%s/simrate_%s.rds", OUT, a))$obs)))
sim$arm <- factor(sim$arm, levels = arms); obs$arm <- factor(obs$arm, levels = arms)
lab <- data.frame(arm = factor(arms, levels = arms), txt = sprintf("D = %.3f, p %s", pkg$KS_D[2:3],
                  ifelse(pkg$KS_p[2:3] < 0.001, "< 0.001", sprintf("= %.2f", pkg$KS_p[2:3]))))
p <- ggplot() +
  geom_density(data = sim, aes(rate, colour = "Predicted", fill = "Predicted"), linewidth = 0.6, alpha = 0.18, adjust = 1.1) +
  geom_density(data = obs, aes(rate, colour = "Observed", fill = "Observed"), linewidth = 0.7, alpha = 0.12, adjust = 1.1, linetype = "22") +
  geom_text(data = lab, aes(x = 62, y = Inf, label = txt), vjust = 1.9, hjust = 1, size = 3.1, colour = "grey20") +
  scale_colour_manual(values = c(Observed = "grey15", Predicted = "#2a78d6"), name = NULL) +
  scale_fill_manual(values = c(Observed = "grey55", Predicted = "#2a78d6"), name = NULL) +
  scale_x_continuous(limits = c(0, 65), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression("Migration rate (km d"^-1*")"), y = "Density") +
  facet_wrap(~arm, ncol = 1) + theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10, hjust = 0), panel.border = element_rect(colour = "grey25", linewidth = 0.4),
        legend.position = "top", legend.key = element_blank(), plot.margin = margin(6, 10, 4, 6))
ggsave(file.path(OUT, "FigS13-1_calibration.jpeg"), p, width = 6.5, height = 4.3, dpi = 400)
cat("S13.1:", lab$txt, "\n")

## ---- S14.1: the by-year tornado, reading the package Table 4 ----
s14 <- readLines("Supplement_Figures/S14_tornado_byyear.R")
cat(tail(run_patched(s14), 45), sep = "\n")

## ---- S10.1: mean-velocity simulation (script 9) and the Figure 4 code (script 11) ----
s9 <- readLines("9_Emm_Eh_sims.R")
a <- grep("^scenarios <- list\\(", s9); b <- a; while(!grepl("\\)\\)\\s*$", s9[b])) b <- b + 1
s9 <- c(s9[seq_len(a - 1)], 'scenarios <- list(Figure = list(landscape = "PcCoastal_VEL_MEAN", duty = MovingFraction$fitted_mean[1]))', s9[(b + 1):length(s9)])
s9 <- fixed(s9, '"Output/Sims/MH_Figure.csv"', sprintf('"%s/MH_Figure_VELMEAN.csv"', OUT))
run_patched(s9)
s11 <- readLines("11_Figures_4_5.R")
s11 <- fixed(s11, 'read.csv("Output/Sims/MH_Figure.csv")', sprintf('read.csv("%s/MH_Figure_VELMEAN.csv")', OUT))
s11 <- fixed(s11, '"Output/Figure4_EnergyUsebyStage.jpeg"', sprintf('"%s/FigS10-1_budget_meanvelocity.jpeg"', OUT))
s11 <- fixed(s11, '"Output/Figure5_Entryday_SpnSucc_byTempUc.jpeg"', 'tempfile(fileext = ".jpeg")')
s11 <- c(s11, sprintf('saveRDS(MH_all, "%s/MH_VELMEAN_budget.rds")', OUT))
run_patched(s11)
vm <- readRDS(file.path(OUT, "MH_VELMEAN_budget.rds"))
vm$top <- with(vm, E_D + E_Mbd + E_Mr + E_H + E_G + E_S)
for(y in c(2022, 2023)){ d <- subset(vm, Year == y & SpawnSiteSuccess == 1)
  w <- subset(d, EntryDate_J >= 209 & EntryDate_J <= 337)
  cat(sprintf("S10.1 %d: max all %.1f | max plotted %.1f | min plotted %.2f on DOY %d | plotted dates <= 9.5: %d of %d | <= 10.3: %d | max success plotted %.1f%% | halted entries plotted %d (last DOY %s) | last entry reaching %d\n",
      y, max(d$top), max(w$top), min(w$top), w$EntryDate_J[which.min(w$top)], sum(w$top <= 9.5), nrow(w), sum(w$top <= 10.3),
      max(w$SpnSuccess_pct), sum(w$Days_halted > 0), max(c(NA, w$EntryDate_J[w$Days_halted > 0]), na.rm = TRUE), max(d$EntryDate_J))) }

## ---- S15.1: EMr split into maintenance while halted, aerobic and anaerobic ----
MF <- readRDS("Output/MovingFraction_levels.rds"); DUTY <- MF$fitted_mean[1]
ID <- read.csv("Output/Landscapes/data_ID_Sacramento.csv")
Spawning_goal <- NULL
s9txt <- readLines("9_Emm_Eh_sims.R")
goal_line <- s9txt[grep("^goal_row", s9txt)]; start_line <- s9txt[grep("^start_row", s9txt)]
fig <- read.csv("Output/Sims/MH_Figure.csv")
OS <- readRDS("Output/OptSwim/OptSwim_tempvel_5.1kg_PcCoastal.rds")
key <- paste(OS$temp, OS$Uc)
set.seed(1)  # as in script 9, where the goal row comes from sampled redd locations
g <- fig$RKM[1]
trace_one <- function(L, entry, goal){
  rown <- 1; coln <- entry; halt <- 0; aer <- 0; ana <- 0
  while(rown < goal & trunc(coln) <= ncol(L$df_t)){
    r <- trunc(rown); cc <- trunc(coln); Um <- as.numeric(L$df_Um[r, cc]); Pm <- as.numeric(L$df_Pm[r, cc]); CO <- as.numeric(L$df_COST[r, cc])
    kmd <- Um*86400/1000*DUTY; md <- ifelse(kmd >= 2, 2, kmd); mv <- md/2
    if(mv == 0){ dt <- 1; halt <- halt + Pm*24*j }else{
      dt <- ifelse((1/kmd*md) >= 1, 1, (1/kmd*md))
      o <- match(paste(L$df_t[r, cc], L$df_v[r, cc]), key); Ps <- as.numeric(OS$Ps[o]); TP <- as.numeric(OS$TravProp[o])
      #df_COST is per elapsed time (script 7): swimming for TP of the time, the rest is recovery
      aer <- aer + (DUTY*TP*(Pm + Ps) + (1 - DUTY)*Pm)*dt*24*j
      ana <- ana + DUTY*(CO - TP*(Pm + Ps))*dt*24*j }
    rown <- rown + mv; coln <- coln + dt
    if(coln >= ncol(L$df_t)) break }
  c(halt = halt, aer = aer, ana = ana) }
goal <- which.min(abs(ID$RKM - 457.3))
res <- bind_rows(lapply(c(2022, 2023), function(y){ L <- readRDS(sprintf("Output/Landscapes/Landscape_PcCoastal_VEL_10pct_%d.rds", y))
  x <- t(sapply(180:365, function(e) trace_one(L, e, goal)))
  data.frame(year = y, entry = 180:365, x) }))
res$total <- with(res, halt + aer + ana)
cmp <- merge(res, fig[, c("Year", "EntryDate_J", "COST")], by.x = c("year", "entry"), by.y = c("Year", "EntryDate_J"))
cat(sprintf("S15.1 check against script 9: max |difference| = %.4f MJ/kg\n", max(abs(cmp$total - cmp$COST*j))))
write.csv(res, file.path(OUT, "FigS15-1_costsplit.csv"), row.names = FALSE)
.jd <- function(d) as.Date(d - 1, origin = "2025-01-01")
lv <- c("Anaerobic overhead", "Aerobic swimming", "Maintenance while unable to advance")
long <- res %>% pivot_longer(c(halt, aer, ana), names_to = "k", values_to = "MJ") %>%
  mutate(component = factor(recode(k, halt = lv[3], aer = lv[2], ana = lv[1]), levels = lv), pdate = .jd(entry), Year = factor(year))
win <- data.frame(Year = factor(c(2022, 2023)), x1 = .jd(c(258, 234)), x2 = .jd(c(329, 292)))
ymax <- ceiling(max(res$total)/2)*2 + 2
p <- ggplot(long, aes(pdate, MJ, fill = component)) +
  geom_rect(data = win, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf), fill = "grey60", alpha = .16) +
  geom_area(colour = "white", linewidth = .1) +
  geom_text(data = win, inherit.aes = FALSE, aes(x = x1 + (x2 - x1)/2, y = ymax - 1.1, label = "observed range of\nentry dates"), size = 2.8, colour = "grey25", lineheight = .95) +
  facet_wrap(~Year) +
  scale_fill_manual(values = c("#D55E00", "#009E73", "#0072B2"), name = NULL) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month", expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, 2), expand = c(0, 0)) +
  labs(x = "Date of passing Rio Vista (RKM 0)", y = expression(paste("Mainstem migration cost, ", E[Mr], " (MJ ", kg^-1, ")"))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = .3), panel.border = element_rect(colour = "grey25", linewidth = .4),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),
        legend.key.size = unit(.9, "lines"), legend.text = element_text(size = 9))
ggsave(file.path(OUT, "FigS15-1_EMr_composition.jpeg"), p, width = 8.6, height = 4.4, dpi = 300)
for(y in c(2022, 2023)){ d <- subset(res, year == y); w <- subset(d, entry >= c(258, 234)[y - 2021] & entry <= c(329, 292)[y - 2021])
  cat(sprintf("S15.1 %d: max total %.1f on DOY %d | halted max %.1f, last halted entry DOY %d | in-window total %.2f-%.2f | in-window anaerobic share max %.0f%% | in-window halted %.2f | max anaerobic share all %.0f%%\n",
      y, max(d$total), d$entry[which.max(d$total)], max(d$halt), max(d$entry[d$halt > 0]), min(w$total), max(w$total), 100*max(w$ana/w$total), max(w$halt), 100*max(d$ana/d$total))) }
