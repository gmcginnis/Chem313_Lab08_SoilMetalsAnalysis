## Week 03 Data Analysis for AA data
## Gillian McGinnis
## Created 06 November 2020
## Updated 06 November 2020

source("R/week02_fullTidy.R")
rm(list = setdiff(ls(), c("AA_merged", "dataKey_tidy")))

dataAA <- AA_merged %>%
  mutate(xRsd = as.numeric(as.character(case_when(xRsd != "HIGH" ~ xRsd))))
dataKey <- dataKey_tidy

rm(list = setdiff(ls(), c("dataAA", "dataKey")))

sampleSites <- unique(filter(dataAA, site != "MB", site != "")$site)

sample_key <- dataKey %>%
  clean_names()

data_aa <- dataAA %>%
  clean_names() %>%
  mutate(rsd = x_rsd)


# Initiate loop, filter out cal data
calAA <- NULL

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
  
equation <- tibble(metal = slope, slopeStd, intercept, interceptStd)
calAA <- rbind(calAA, equation)

remove(equation, cal, slope, slopeStd, intercept, interceptStd, w, model, uniqueMetal)