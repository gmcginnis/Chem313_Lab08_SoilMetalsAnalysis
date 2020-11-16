## Week 03 Data Analysis for AA data
## Gillian McGinnis
## Created 11 November 2020
## Updated 15 November 2020

source("R/week03_AAanalysis.R")
##units: ppm
library(broom)

stats_AA <- sample_data %>%
  group_by(site) %>%
  summarize(mean_conc = mean(conc_blanked),
            sd_conc = sd(conc_blanked),
            n = n()) %>%
  mutate(se = qnorm(0.975)*sd_conc/sqrt(n),
         lower_ci = mean_conc - se,
         upper_ci = mean_conc + se)

#remove(sample_data)

## Stat tests:
sample_sites <- unique(sample_data$site)
metals_analyzed <- unique(sample_data$metal)

filtered_AA <- sample_data %>%
  filter(site != "QC")

## ANOVA:  
anova <- aov(conc_blanked ~ site, data = filtered_AA)
anova_df <- anova %>%
  tidy() %>%
  as.data.frame() %>%
  mutate(label = case_when(
    p.value < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
    p.value >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
  ))


## Tukey's Post-Hoc HSD:
tukey_table <- TukeyHSD(anova)
plot(tukey_table, las = 1)
#tukey_df <- rbind(tukey_df, as.data.frame(tukey_table$site))
tukey_df <- as.data.frame(tukey_table$site) %>%
  rownames_to_column() %>%
  mutate(label = case_when(
    `p adj` < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
    `p adj` >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
  )) %>%
  rename(pair = rowname)

ggplot(tukey_df, aes(color = label))+
  geom_hline(yintercept=0, lty="11", color="grey30") +
  geom_errorbar(aes(pair, ymin=lwr, ymax=upr), width=0.2) +
  geom_point(aes(pair, diff)) +
  labs(color="",
       x = "Location pairing",
       y = "Difference")+
  scale_x_discrete()+
  theme_few()+
  coord_flip()


## T-tests
pairwise.t.test(filtered_AA$conc_blanked, filtered_AA$site)

t.test(conc_blanked ~ site, data=subset(filtered_AA, site %in% c("C", "F")))
