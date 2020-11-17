## Week 03 Data Analysis for ICPMS data
## Gillian McGinnis
## Created 11 November 2020
## Updated 15 November 2020

source("R/week03_ICPMSanalysis.R")
library(broom)
#Units: ppb
stats_ICPMS <- sample_data %>%
  group_by(metal, site) %>%
  summarize(mean_conc = mean(conc_blanked),
            sd_conc = sd(conc_blanked),
            n = n()) %>%
  mutate(se = qnorm(0.975)*sd_conc/sqrt(n),
         lower_ci = mean_conc - se,
         upper_ci = mean_conc + se) %>%
  mutate(se_ppm = se/1000,
         lower_ci_ppm = lower_ci/1000,
         upper_ci_ppm = upper_ci/1000)

ICPMS_errors <- sample_data %>%
  group_by(metal, site) %>%
  summarize(mean_conc = mean(conc_blanked),
            sd_conc = sd(conc_blanked),
            mean_conc_error = mean(conc_blanked_error)) %>%
  mutate(mean_conc = mean_conc/1000,
         mean_conc_error = mean_conc_error/1000)

ICPMS_error_prop <- sample_data %>%
  select(!c(conc_unblanked, conc_unblanked_error)) %>%
  group_by(site, metal) %>%
  mutate(numerator = (conc_blanked-conc_blanked_error)^2) %>%
  summarize(num_sum = sum(numerator)) %>%
  mutate(n = case_when(
    site == "A" ~ 4,
    site == "B" ~ 5,
    site == "C" ~ 5,
    site == "D" ~ 3,
    site == "E" ~ 2,
    site == "F" ~ 3,
    site == "QC" ~ 12
  )) %>%
  group_by(site,metal) %>%
  summarize(error_prop = sqrt(num_sum/(n-1))) %>%
  mutate(error_prop_ppm = error_prop/1000)
  #mutate(error_prop = sqrt((conc_blanked-conc_blanked_error)^2)/(n-1))
  # summarize(mean_conc = mean(conc_blanked),
  #           n_new = n(),
  #           error_prop = sqrt((conc_blanked-conc_blanked_error)^2)/(n_new-1))

ICPMS_minimal_error_prop <- sample_data %>%
  mutate(metal_short= case_when(
    metal == "As75" ~ "As",
    metal == "Cd111" ~ "Cd",
    metal == "Cd114" ~ "Cd",
    metal == "Cr52" ~ "Cr",
    metal == "Cr53" ~ "Cr",
    metal == "Pb208" ~ "Pb"
  )) %>%
  group_by(metal_short, site) %>%
  select(!c(conc_unblanked, conc_unblanked_error)) %>%
  group_by(site, metal_short) %>%
  mutate(numerator = (conc_blanked-conc_blanked_error)^2) %>%
  summarize(num_sum = sum(numerator)) %>%
  mutate(n = case_when(
    site == "A" ~ 4,
    site == "B" ~ 5,
    site == "C" ~ 5,
    site == "D" ~ 3,
    site == "E" ~ 2,
    site == "F" ~ 3,
    site == "QC" ~ 12
  )) %>%
  group_by(site,metal_short) %>%
  summarize(error_prop = sqrt(num_sum/(n-1))) %>%
  mutate(error_prop_ppm = error_prop/1000)

minimal_stats_ICPMS <- sample_data %>%
  mutate(metal_short= case_when(
    metal == "As75" ~ "As",
    metal == "Cd111" ~ "Cd",
    metal == "Cd114" ~ "Cd",
    metal == "Cr52" ~ "Cr",
    metal == "Cr53" ~ "Cr",
    metal == "Pb208" ~ "Pb"
  )) %>%
  group_by(metal_short, site) %>%
  summarize(mean_conc = mean(conc_blanked),
            sd_conc = sd(conc_blanked),
            n = n()) %>%
  mutate(se = qnorm(0.975)*sd_conc/sqrt(n),
         lower_ci = mean_conc - se,
         upper_ci = mean_conc + se) %>%
  mutate(mean_ppm = mean_conc/1000,
         sd_ppm = sd_conc/1000,
         se_ppm = se/1000,
         lower_ci_ppm = lower_ci/1000,
         upper_ci_ppm = upper_ci/1000) %>%
  select(site, metal_short, mean_ppm, sd_ppm, se_ppm, lower_ci_ppm, upper_ci_ppm)
  # mutate(se = sd_conc/sqrt(n),
  #        lower_ci = mean_conc - qt(1 - (0.05/2), n - 1) * se,
  #        upper_ci = mean_conc + qt(1 - (0.05/2), n - 1) * se)
#remove(sample_data)

all_ICPMS_stats <- full_join(minimal_stats_ICPMS, ICPMS_minimal_error_prop) %>%
  select(site, metal_short, mean_ppm, sd_ppm, error_prop_ppm)

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
  mutate(mass_fraction = mean_conc/1000,
         mass_fraction_sd = sd_conc/1000) %>%
  mutate(metal_short= case_when(
    metal == "As75" ~ "As",
    metal == "Cd111" ~ "Cd",
    metal == "Cd114" ~ "Cd",
    metal == "Cr52" ~ "Cr",
    metal == "Cr53" ~ "Cr",
    metal == "Pb208" ~ "Pb"
  ))

joining_given <- given_qc %>%
  rename(metal_short = metal,
         mass_frac_given = mass_frac,
         mass_frac_sd_given = mass_frac_sd) %>%
  select(!mass_frac_sd_given)
  # rename(metal = metal,
  #        mass_frac = mass_frac,
  #        mass_frac_sd = mass_frac_sd)

joining_ICPMS <- stats_qc %>%
  rename(mass_frac_icpms = mass_fraction,
         mass_frac_sd_icpms = mass_fraction_sd) %>%
  select(c(metal, mass_frac_icpms, metal_short))
  # rename(mass_frac = mass_fraction,
  #        mass_frac_sd = mass_fraction_sd) %>%
  # select(c(metal, mass_frac, mass_frac_sd))

joined_qc <- full_join(joining_given, joining_ICPMS) %>%
  drop_na() %>%
  mutate(per_recovery = (mass_frac_icpms/mass_frac_given)*100)
## Df with percent recovery ^

# t_df <- NULL
# #t_test_function <- function(metals){
#   tidied <- tidy(t.test(mass_frac ~ metal, data=subset(joined_qc, metal %in% c(metals))), ) %>%
#     as.data.frame() %>%
#     mutate(pair = paste(metals, collapse=""))
#   
#   t_df <<- rbind(t_df, tidied)
# #}
# t_test_function(c("Cd", "Cd111"))
# metals <- c("Cd", "Cd111")
# subset_metals <- subset(joined_qc, metal %in% c(metals))


## Stat tests
sample_sites <- unique(sample_data$site)
metals_analyzed <- unique(sample_data$metal)

anova_df <- NULL

aov_test <- function(unique_metal){
  filtered_ICPMS <- sample_data %>%
    filter(site != "QC") %>%
    filter(metal == unique_metal)
  
  anova <- aov(conc_blanked ~ site, data = filtered_ICPMS) %>%
    tidy()
  anova <- as.data.frame(anova) %>%
    mutate(metal = unique_metal)
  anova_df <<-rbind(anova_df, anova)
  #return(anova_df)
}

# I know I could write a function for this but it's being finicky
aov_test("As75")
aov_test("Cd111")
aov_test("Cd114")
aov_test("Cr52")
aov_test("Cr53")
aov_test("Pb208")

anova_df <- anova_df %>%
  mutate(label = case_when(
    p.value < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
    p.value >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
  ))
# Surprise! Difference is not statistically significant for all metals.

## Tukey's Post-Hoc HSD
# I made the mistake of doing this BEFORE checking my ANOVA results... they are all non-sig
tukey_df <- NULL

tukey_test <- function(unique_metal){
  filtered_ICPMS <- sample_data %>%
    filter(site != "QC") %>%
    filter(metal == unique_metal)
  
  anova <- aov(conc_blanked ~ site, data = filtered_ICPMS)
  tukey_table <- TukeyHSD(anova)
  plot(tukey_table, las = 1, sub = paste(unique_metal))
  #tukey_df <- rbind(tukey_df, as.data.frame(tukey_table$site))
  tukey_frame <- as.data.frame(tukey_table$site) %>%
    rownames_to_column() %>%
    mutate(metal = unique_metal)
  
  tukey_df <<- rbind(tukey_df, tukey_frame)
  #print(tukey_df)
  #return(tukey_df)
}

tukey_df <- NULL
# I know I could write a function to run all the metals but it's being finicky
tukey_test("As75")
tukey_test("Cd111")
tukey_test("Cd114")
tukey_test("Cr52")
tukey_test("Cr53")
tukey_test("Pb208")

tukey_df <- tukey_df %>%
  mutate(label = case_when(
    `p adj` < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
    `p adj` >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
    )) %>%
  rename(pair = rowname)

ggplot(tukey_df, aes(color = label))+
  facet_wrap(~metal)+
  geom_hline(yintercept=0, lty="11", color="grey30") +
  geom_errorbar(aes(pair, ymin=lwr, ymax=upr), width=0.2) +
  geom_point(aes(pair, diff)) +
  labs(color="",
       x = "Location pairing",
       y = "Difference")+
  scale_x_discrete()+
  theme_few()+
  coord_flip()
