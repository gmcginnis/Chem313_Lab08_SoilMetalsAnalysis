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
         upper_ci = mean_conc + se)
  # mutate(se = sd_conc/sqrt(n),
  #        lower_ci = mean_conc - qt(1 - (0.05/2), n - 1) * se,
  #        upper_ci = mean_conc + qt(1 - (0.05/2), n - 1) * se)
#remove(sample_data)

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
