## Week 03 Data Analysis for ICPMS data
## Gillian McGinnis
## Created 11 November 2020
## Updated 11 November 2020

source("R/week03_ICPMSanalysis.R")
#Units: ppb
stats_ICPMS <- sample_data %>%
  group_by(metal, site) %>%
  summarize(mean_conc = mean(conc_blanked),
            sd_conc = sd(conc_blanked),
            n = n()) %>%
  mutate(se = qnorm(0.975)*sd_conc/sqrt(n),
         lower_ci = mean_conc - se,
         upper_ci = mean_conc + se)
  # mutate(se = sd_conc/sqrt(n),
  #        lower_ci = mean_conc - qt(1 - (0.05/2), n - 1) * se,
  #        upper_ci = mean_conc + qt(1 - (0.05/2), n - 1) * se)
remove(sample_data)

stats_ICPMS_Cr <- stats_ICPMS %>%
  filter(metal == "Cr53",
         site != "QC")

#Units: mass fraction (mg/kg)
given_qc <- data.frame(
  metal = c("Cd", "Cr", "Pb"),
  mass_frac = c(2.94, 121.9, 150),
  mass_frac_sd = c(0.29, 3.8, 17)
)

stats_qc <- stats_ICPMS %>%
  filter(site == "QC") %>%
  mutate(mass_fraction = mean_conc/1000)
