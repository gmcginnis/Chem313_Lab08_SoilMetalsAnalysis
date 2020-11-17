## Full tidying code. Copied from week02_ICPMStidy.Rmd and week02_AAtidy.R
# Used for week03 data analysis
# Doesn't include writing new csv files.

library(tidyverse)
library(readr)
library(janitor)


## ICPMS
# Importing necessary data
ICPMS_imported <- read.csv("data/ICPMS_Data.csv", skip = 1, na = "N/A")
sample_key <- read.csv("data/Sample_Key.csv", skip = 0)

RSD_data <- ICPMS_imported %>%
  # picking relevant columns + renaming them to appropraite isotopes
  select(Cr52 = CPS.RSD,
         Cr53 = CPS.RSD.1,
         As75 = CPS.RSD.2,
         Cd111 = CPS.RSD.3,
         Cd114 = CPS.RSD.4,
         Pb208 = CPS.RSD.5,
         Ge_RSD = CPS.RSD.7,
         Sample.Key) %>%
  # making data longer
  pivot_longer(1:6,
               names_to = "metal",
               values_to = "RSD")

ICPMS_tidy <- ICPMS_imported %>%
  select(Cr52 = CPS,
         Cr53 = CPS.1,
         As75 = CPS.2,
         Cd111 = CPS.3,
         Cd114 = CPS.4,
         Pb208 = CPS.5,
         Ge72 = CPS.7,
         Sample.Key) %>%
  pivot_longer(1:6,
               names_to = "metal",
               values_to = "CPS") %>%
  mutate(RSD = RSD_data$RSD/RSD_data$Ge_RSD,
         CPS = CPS/Ge72) %>%
  select(-Ge72)

# Confirming the RSD data matches
all(RSD_data$Sample.Key == ICPMS_tidy$Sample.Key,
    RSD_data$metal == ICPMS_tidy$metal)

ICPMS_merged <- merge(ICPMS_tidy, sample_key)



## AA

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
ggplot(AA_merged, aes(x = concentration, y = meanAbs))+
  geom_point()+
  stat_smooth(method = "lm")