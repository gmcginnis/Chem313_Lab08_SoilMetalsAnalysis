## Week 02 AA Tidy
## Gillian McGinnis
## Created 20 October 2020
## Updated 20 October 2020

## Loading necessary libraries
library(tidyverse)
library(readr)
library(janitor)

## Importing dataframes
dataAA <- read.csv("data/AA_Data.csv", skip = 4)
dataKey <- read.csv("data/Sample_Key.csv")

## Tidying data
dataKey_tidy <- dataKey %>%
  clean_names("small_camel") %>%
  mutate(sampleKey = as.character(sampleKey))
  
dataAA_tidy <- dataAA %>%
  clean_names("small_camel")

dataMerged <- full_join(dataKey_tidy, dataAA_tidy)

AA_merged <- dataMerged %>%
  drop_na(meanAbs)

## Quick plot to check
ggplot(AA_merged, aes(x = meanAbs, y = concentration))+
  geom_point()+
  stat_smooth(method = "lm")


saveRDS(AA_merged, file = "data/AA_merged.rds")
write.csv(AA_merged, file = "data/AA_merged.csv")

## Test code below
# dataAA_tidy <- dataAA %>%
#   clean_names("small_camel") %>%
#   mutate(tag = case_when(
#     sampleKey ==  "Sample Blank" | sampleKey == "check10" ~ "Standard",
#     sampleKey != "Sample Blank" & sampleKey != "check10" ~ "Sample"))
# %>%mutate(altKey = case_when(
#   sampleKey ==  "Sample Blank" | sampleKey == "check10" ~ 00,
#   sampleKey != "Sample Blank" & sampleKey != "check10" ~ sampleKey))
  

# dataMerged <- merge(dataAA, dataKey)
# dataMerged <- dataAA %>%
#   merge(dataKey) %>%
#   clean_names(case = "small_camel")
#   select(-c(Mean.Abs., X.RSD))