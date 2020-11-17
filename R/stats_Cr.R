## Week 03 Data Analysis for Cr data, comparing ICPMS and AA
## Gillian McGinnis
## Created 15 November 2020
## Updated 16 November 2020

#use Cr53
source("R/week03_ICPMSanalysis.R")
sample_data_ICPMS <- sample_data
remove(sample_data)
save(data = sample_data_ICPMS, file = "data/sample_data_ICPMS.Rdata")

source("R/week03_AAanalysis.R")
sample_data_AA <- sample_data
remove(sample_data)
save(data = sample_data_AA, file = "data/sample_data_AA.Rdata")

rm(list = ls())

load("data/sample_data_ICPMS.Rdata")
load("data/sample_data_AA.Rdata")

cr_icpms <- sample_data_ICPMS %>%
  filter(metal == "Cr53") %>%
  select(site, conc_blanked, conc_blanked_error) %>%
  mutate(conc_blanked = conc_blanked/1000,
         conc_blanked_error = conc_blanked_error/1000) %>%
  mutate(instrument = "ICPMS")

cr_aa <- sample_data_AA %>%
  select(site, conc_blanked, conc_blanked_error) %>%
  mutate(instrument = "AA")

cr_all <- full_join(cr_icpms, cr_aa)

sample_sites <- unique(cr_all$site)
#Input: location
# Output: df of t.test results
t_df <- NULL
t_test_function <- function(site_input){
  tidied <- tidy(t.test(conc_blanked ~ instrument, data = subset(cr_all, site == site_input))) %>%
    as.data.frame() %>%
    mutate(site = site_input)
  
  t_df <<- rbind(t_df, tidied)
}

run_sites <- function(Function){
  value <- NULL
  for(sites in sample_sites){
    site_value <- Function(sites)
    value <- rbind(site_value, value)
  }
  return(value)
}

run_sites(t_test_function)

t_df <- t_df %>%
  mutate(label = case_when(
    p.value < 0.05 ~ "p < 0.05", # Reject null hypothesiss; diff is significant
    p.value >= 0.05 ~ "Non-Sig" # Fail to reject null hyp; diff is not significant
  ))
