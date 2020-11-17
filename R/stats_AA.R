## Week 03 Data Analysis for AA data
## Gillian McGinnis
## Created 11 November 2020
## Updated 16 November 2020

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

# ## unnecessary error prop
# AA_error_means <- sample_data %>%
#   group_by(site) %>%
#   summarize(mean = mean(conc_blanked))
# 
# stats_AA_error <- sample_data %>%
#   full_join(AA_error_means) %>%
#   group_by(site) %>%
#   mutate(numerator = (conc_blanked-mean)^2) %>%
#   #drop_na(numerator) %>%
#   summarize(num_sum = sum(numerator)) %>%
#   mutate(n = case_when(
#     site == "A" ~ 4,
#     site == "B" ~ 5,
#     site == "C" ~ 5,
#     site == "D" ~ 3,
#     site == "E" ~ 2,
#     site == "F" ~ 3,
#     site == "QC" ~ 12
#   )) %>%
#   group_by(site) %>%
#   summarize(error_prop = sqrt(num_sum/(n-1)))
# 
# # stats_AA_error <- sample_data %>%
# #   group_by(site) %>%
# #   mutate(numerator = (conc_blanked-conc_blanked_error)^2) %>%
# #   drop_na(numerator) %>%
# #   summarize(num_sum = sum(numerator)) %>%
# #   mutate(n = case_when(
# #     site == "A" ~ 4,
# #     site == "B" ~ 5,
# #     site == "C" ~ 5,
# #     site == "D" ~ 3,
# #     site == "E" ~ 2,
# #     site == "F" ~ 3,
# #     site == "QC" ~ 12
# #   )) %>%
# #   group_by(site) %>%
# #   summarize(error_prop = sqrt(num_sum/(n-1)))
# 
# all_AA_stats <- full_join(stats_AA, stats_AA_error) %>%
#   select(site, metal_short, mean_ppm, sd_ppm, error_prop_ppm)



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
# sad

## T-tests
#Guess it's time to go old school:

t_df <- NULL
t_test_function <- function(pairs){
  tidied <- tidy(t.test(conc_blanked ~ site, data=subset(filtered_AA, site %in% c(pairs)))) %>%
    as.data.frame() %>%
    mutate(pair = paste(pairs, collapse=""))
  
  t_df <<- rbind(t_df, tidied)
}
#Again I'm sure there's a way to loop this but here we are
t_test_function(c("A", "B"))
t_test_function(c("A", "C"))
t_test_function(c("A", "D"))
t_test_function(c("A", "E"))
t_test_function(c("A", "F"))
t_test_function(c("B", "C"))
t_test_function(c("B", "D"))
t_test_function(c("B", "E"))
t_test_function(c("B", "F"))
t_test_function(c("C", "D"))
t_test_function(c("C", "E"))
t_test_function(c("C", "F"))
t_test_function(c("D", "E"))
t_test_function(c("D", "F"))
t_test_function(c("E", "F"))

t_df <- t_df %>%
  mutate(label = case_when(
    p.value < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
    p.value >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
  ))

#Neat, right??

#t_df$pair <- factor(t_df$pair, levels = c('AB', 'AC', 'AD', 'AE', 'AF', 'BC', 'BD', 'BE', 'BF', 'CD', 'CE', 'CF', 'DE', 'DF', 'EF'))
t_df$pair <- factor(t_df$pair, levels = c('EF','DF','DE' ,'CF' ,'CE' ,'CD' ,'BF' ,'BE' ,'BD' ,'BC' ,'AF' ,'AE' ,'AD' ,'AC' ,'AB'))

t_test_viz <- ggplot(t_df, aes(color = label))+
  geom_hline(yintercept=0, lty="11", color="grey30") +
  geom_errorbar(aes(pair, ymin=conf.low, ymax=conf.high), width=0.2) +
  geom_point(aes(pair, parameter)) +
  labs(color="",
       x = "Location pairing",
       y = "CI")+
  scale_x_discrete()+
  theme_few()+
  coord_flip()

ggsave("t_test_AA.png", plot = t_test_viz, path = "figures/")


paired_t_test <- pairwise.t.test(filtered_AA$conc_blanked, filtered_AA$site) %>%
  tidy() %>%
  as.data.frame() %>%
  mutate(label = case_when(
    p.value < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
    p.value >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
  ))

unpaired_t_test <- filtered_AA %>%
  group_by(site) %>%
  do(tidy(t.test(.$conc_blanked))) %>%
  mutate(label = case_when(
    p.value < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
    p.value >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
  ))


filtered_AA %>%
  group_by(site) %>%
  doo(~t.test(conc_blanked ~ site, data=.)) %>%
  mutate(label = case_when(
    p.value < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
    p.value >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
  ))



# paired_t_test <- filtered_AA %>%
#   do(tidy(t.test(.$conc_blanked ~ .$site))) %>%
#   mutate(label = case_when(
#     p.value < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
#     p.value >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
#   ))

# stat_test <- filtered_AA %>%
#   t_test(conc_blanked ~ site)

## LIMITS:

mean_AA <- mean(filtered_AA$conc_blanked)
sd_AA <- sd(filtered_AA$conc_blanked)

## limit of detection
## dl = 3s/m
lod <- (3*sd_AA)/mean_AA
lod
## limit of quantification
## ql = 10s/m
loq <- (10*sd_AA)/mean_AA
loq