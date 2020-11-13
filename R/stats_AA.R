## Week 03 Data Analysis for ICPMS data
## Gillian McGinnis
## Created 11 November 2020
## Updated 11 November 2020

source("R/week03_AAanalysis.R")
##units: ppm

stats_AA <- sample_data %>%
  group_by(site) %>%
  summarize(mean_conc = mean(conc_blanked),
            sd_conc = sd(conc_blanked),
            n = n()) %>%
  mutate(se = qnorm(0.975)*sd_conc/sqrt(n),
         lower_ci = mean_conc - se,
         upper_ci = mean_conc + se)

remove(sample_data)