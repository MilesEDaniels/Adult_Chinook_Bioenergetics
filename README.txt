Code and data for:
"Modelling temperature and velocity effects on migration energetics and spawning success of adult
Chinook salmon in a regulated river" (Conservation Physiology)
Alyssa FitzGerald, Eric Danner, Benjamin Atencio, Nate Mantua, Cyril Michel, Christopher C. Caudill,
Jeremy Notch, Miles E. Daniels
Contact: Miles E. Daniels, miles.daniels@ucsc.edu

The following describes the scripts required to reproduce the tables and figures of the main text,
with a brief description for each. The energy budget model is adapted from FitzGerald & Martin (2022,
Conservation Physiology) and uses the bioenergetics model (Model 3) of Martin et al. (2015, Ecology Letters).
Figures 1 (conceptual diagram) and 2 (map) were not made in R.

HOW TO RUN
Set the R working directory to this folder and run the scripts in order, 1 to 12.
Each script reads from Data/ and from the Output/ files written by earlier scripts, and writes to Output/.
Total run time is about 4 minutes on a laptop (script 10 takes about 3 minutes).
Tested with R 4.4.2 and packages ggplot2 3.5.1, dplyr 1.1.4, ggtext 0.1.2, multcompView 0.1.10,
MASS 7.3-61, truncnorm 1.0.9, cowplot 1.1.3, colorRamps 2.3.4, scales 1.3.0. Figures were saved with
the ragg graphics device (ragg 1.5.0), which ggsave() uses when it is installed.

SCRIPTS
1. 1_Fatmeter_EnergyDensity.R: Convert fatmeter scores to somatic energy density; summary statistics by
   migration stage; ANOVA and Tukey HSD tests between stages and between years within stages;
   Table 1, Figure 3; EI and ED in Table 2.
2. 2_FL_Mass.R: Model mass from fork length (Klamath River adults; Hearsey 2011); Supplement 1; used in script 5.
3. 3_Eg_allocation.R: Model the percent of initial somatic energy that females allocate to gonad
   formation (EG); Table 2, Supplement 4.
4. 4_Es_parameterization.R: Model the energy density used for spawning (ES); must run script 1 first;
   Table 2, Supplement 3.
5. 5_Embd_BayDelta.R: Estimate energy used by tagged fish migrating through the bay and delta (EMbd),
   assuming tidal surfing or active swimming, and fit its distribution; summarise the date fish passed
   Rio Vista, used to weight entry timing; must run script 2 first; Table 2, Supplement 7.
6. 6_OptimalMigSwimSpeed.R: Predict the swim speed that minimises cost of transport, and the resulting
   migration speed, for each temperature and water velocity, under each formulation of critical power
   (maximum aerobic capacity): the base model, and the aerobic scope of the coastal fall-run population
   (Van Wert et al. 2023) with its upper pejus temperature at the mean and +/- 1 SD; Supplement 5.
7. 7_EnergyLandscapes.R: Create temperature, velocity, migration speed, maintenance metabolism and cost
   landscapes (RKM x day matrices) for the mainstem Sacramento River in 2012, 2022 and 2023, for each
   critical power formulation and velocity metric (10th percentile or mean); must run script 6 first.
8. 8_MovingFraction_Calibration.R: Calibrate the moving fraction (the fraction of each time step a fish
   migrates): its shape on fish tracked in 2012 and its mean on fish tracked in 2022-2023; compare
   observed and predicted migration rates (RKM 20 to 150) with KS tests; Results, Table 3;
   must run script 7 first.
9. 9_Emm_Eh_sims.R: Simulate energy used during mainstem migration (EMr) and holding (EH) for one fish
   entering freshwater on every day of year 180-365 in 2022 and 2023, for each level of the moving
   fraction and each sensitivity scenario; must run scripts 7 and 8 first.
10. 10_AdultSurvival_Sensitivity.R: Use the energy budget to estimate the likelihood that each fish
    spawns successfully, weighted by observed entry timing; baseline spawning success and the
    one-at-a-time sensitivity analysis; Results, Table 4; must run script 9 first.
11. 11_Figures_4_5.R: Plot energy use by stage (Figure 4) and predicted spawning success by date of
    freshwater entry, coloured by temperature and velocity exposure (Figure 5); must run script 9 first.
12. 12_Figure_6.R: Plot the one-at-a-time sensitivity analysis (Figure 6); must run script 10 first.

SUPPLEMENT FIGURES (optional, run after scripts 1-12)
Supplement_Figures/make_supplement_figures.R: Regenerate Figures S10.1, S13.1, S14.1 and S15.1. Run with
   source("Supplement_Figures/make_supplement_figures.R") from this folder; also needs the tidyr package.
   It writes only to Supplement_Figures/ and takes under a minute. S14_tornado_byyear.R in that folder is
   called by it.

DATA (Data/)
Fatmeter_EnergyDensity_FallRun.csv: one row per fall-run Chinook salmon sampled with the fatmeter,
   2021-2024 (n = 3,954). Study_Year; Location_abrev (sampling site; STAN = Stanislaus River);
   Energy_stage (migration stage); Sex_MF (sex recorded in the field); FL_mm (fork length, mm);
   S1, S2 (the two fatmeter scans).
Klamath_Hearsey2011_FL_Mass_Gonad.csv: Klamath River adult Chinook salmon from Hearsey (2011),
   Appendices F and G. Run, Sex, Stage (Entry or Spawn-ready), FL_mm, Mass_kg, Gonad_mass_g.
Telemetry_BayDelta_2022_2023.csv: acoustically tagged fish that reached Freeport in 2022 (n = 11) and
   2023 (n = 56). Fork_length_cm; GG_date_UTC (first detection past the Golden Gate); RV_date_UTC (first
   detection near, but not past, Rio Vista); GG_dist_km and RV_dist_km (distance from the Golden Gate
   at those detections); Temp_mean_C (mean temperature experienced between them).
Telemetry_RioVista_passage_2022_2023.csv: date each tagged fish was first detected past Rio Vista,
   2022 (n = 16) and 2023 (n = 64). Used for the observed-detections bar in Figure 4.
Telemetry_MigrationRate_RKM20_150.csv: tagged fish detected on the mainstem Sacramento River at both
   RKM 20 and RKM 150 (RKM 0 = Rio Vista), 2012 (n = 18), 2022 (n = 5) and 2023 (n = 46).
   RKM20_date_UTC and RKM150_date_UTC: first time each fish was positioned at or above each RKM, from
   6-hourly interpolated positions. 2012 fish are from Michel et al. (unpublished).
   Rows are in the order of the original analysis; the simulated rates in Table 3 pair fish with random
   draws by position, so changing the row order changes the KS statistics slightly.
River/: daily water temperature (C) and cross-sectional mean and 10th percentile velocity (m/s; files
   ending _10p) of the Sacramento River from Rio Vista to Keswick Dam, from the RAFT model, at 2 km
   resolution, 2012, 2013, 2022 and 2023 (rows = river distance, columns = days).
   River_Distance_Downstream_Keswick_m_2014.csv gives the distance of each row downstream of Keswick Dam (m);
   River_Time_2016.csv gives the dates of the columns (a leap year, used for all years).
Telemetry files are summaries derived from the tracks of individual fish; tag codes are not included.
Fish_ID numbers are arbitrary and do not link fish between files.

FIGURE FILES (Fig_1.pdf ... Fig_6.pdf)
Publication versions of the six main-text figures: each is the figure image placed on a page 6.5 in wide.
Figures 3-6 are the Output/ JPEGs (effective resolution 415-615 dpi); Figures 1 (conceptual diagram) and 2
(map) were not made in R and are included as submitted (300 dpi).

OUTPUT (Output/), written by the scripts
Table1_EnergyDensity_byStage.csv          Table 1 (script 1)
Table3_MigrationRate_calibration.csv      Table 3 (script 8)
Table4_Sensitivity_OAT.csv                Table 4; row "0" is the baseline spawning success (script 10)
Figure3 ... Figure6 .jpeg                 Figures 3-6 (scripts 1, 11, 12)
Eg_allocation.csv, Es_parameter.csv, Embd_and_entry_parameters.csv, FL_Mass_coefficients.csv
                                          Table 2 parameters (scripts 2-5)
MovingFraction_levels.csv                 calibrated moving fraction (script 8)
Survival_bylevel.csv                      spawning success for each level of the moving fraction (script 10)
OptSwim/, Landscapes/, Sims/              intermediate files

NOTES
- Spawning cost (ES): the energy budget samples ES from a normal distribution truncated at zero with
  mean 1.30 and SD 0.84 MJ/kg (scripts 10 and 11), the values estimated by script 4 and reported in Table 2.
- Figure 3: Table 1 and the significance tests include 575 Stanislaus River fish recorded as stage 2 and
  reclassified as stage 3 (spawn-ready) because of their low fatmeter scores. Figure 3 plots stage 3
  without these fish, as published.
- Script 6 corrects an error in the code used during the review, in which the anaerobic cost term was
  evaluated for only the slowest candidate swim speed when choosing the swim speed that minimises cost of
  transport. Results in the article were produced with the corrected code.
- Script 7 corrects a second error in the code used during the review. Model 3 cost (COST) is energy per unit
  time spent swimming, but the simulation charged it per unit elapsed time, which counted recovery twice
  wherever a fish had to stop to recover (TravProp < 1). COST is now multiplied by TravProp, so that energy
  per km equals the cost of transport. Migration timing is unchanged.
