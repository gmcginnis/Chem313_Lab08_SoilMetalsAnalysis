## Week 03 Data Analysis for AA data
## Gillian McGinnis
## Created 06 November 2020
## Updated 06 November 2020

source("R/week02_fullTidy.R")
rm(list = setdiff(ls(), c("AA_merged", "dataKey_tidy")))

dataAA <- AA_merged %>%
  mutate(rsd = case_when(xRsd != "HIGH" ~ xRsd)) %>%
  select(-xRsd) %>%
  mutate(xRsd = as.numeric(as.character(rsd))) %>%
  select(-rsd)

dataAA2 <- AA_merged %>%
  mutate(xRsd = as.numeric(as.character(case_when(xRsd != "HIGH" ~ xRsd))))

dataKey <- dataKey_tidy
rm(list = setdiff(ls(), c("dataAA", "dataKey")))
