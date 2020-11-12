## Week 03 Data Analysis for AA data
## Gillian McGinnis
## Created 06 November 2020
## Updated 11 November 2020

library(tidyverse)
library(ggthemes)
source("R/week02_fullTidy.R")
rm(list = setdiff(ls(), c("AA_merged", "dataKey_tidy")))

## Adjusting dataframes
dataAA <- AA_merged %>%
  mutate(xRsd = as.numeric(as.character(case_when(xRsd != "HIGH" ~ xRsd))))
dataKey <- dataKey_tidy

## Cleaning environment
rm(list = setdiff(ls(), c("dataAA", "dataKey")))

## Creating a list of sample sites by excluding NAs and Method Blanks
sampleSites <- unique(filter(dataAA, site != "MB", site != "")$site)

## Changing names back to snake_case
sample_key <- dataKey %>%
  clean_names()

## Chancing names back to snake_case and renaming a variable to match that of other similar dataframes
data_aa <- dataAA %>%
  clean_names() %>%
  mutate(rsd = x_rsd)


## Code based on ICPMS analysis
# Initiate loop, filter out cal data
calAA <- NULL
## Selecting cal standards
cal <- data_aa %>%
  filter(type == "CalStd" | type == "CalStd2" | type =="CalStd4") %>%
  select(concentration, mean_abs, rsd)
  
  #Step 3.2, perform weighted linear reg, and pull out relevant info from model
w <- 1/ (cal$mean_abs * cal$rsd)^2
model <- lm(cal$mean_abs ~ cal$concentration, weights = w)
#model_unweighted <- lm(cal$mean_abs ~ cal$concentration)
  
slope <- model$coefficients[2]
intercept <- model$coefficients[1]
slopeStd <- summary(model)$coefficients[2,2]
interceptStd <- summary(model)$coefficients[1,2]

plot(cal$mean_abs ~ cal$concentration,
     xlab = paste("Concentration of Cr (ppb)"),
     ylab = "Counts per second")+
  abline(model, col = "red")+
  title(paste("Calibration for Cr"))

## Alternatively, here's my way.
myResults <- summary(lm(mean_abs ~ concentration, weights = w, data = cal))

ggplot(cal, aes(x = concentration, y = mean_abs))+
  geom_point(shape = 1)+
  geom_smooth(method = "lm", se = FALSE, color = "red")+
  theme_few()+
  labs(x = "Concentration of Cr (ppb)", y = "Counts per second", title = "Calibration curve for Cr")
  
## Back to some given stuff:
equation <- tibble(metal = slope, slopeStd, intercept, interceptStd)
calAA <- rbind(calAA, equation)

remove(equation, cal, slope, slopeStd, intercept, interceptStd, w, model, uniqueMetal)


## Step 4.1: creating a function to analyze the samples
