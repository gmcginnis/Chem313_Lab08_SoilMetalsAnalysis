## Week 03 Data Analysis
## Gillian McGinnis
## Created 06 November 2020
## Updated 06 November 2020

source("R/week02_fullTidy.R")
rm(list = setdiff(ls(), c("AA_merged", "ICPMS_merged")))

dataICPMS <- ICPMS_merged %>%
  clean_names("small_camel")
dataAA <- AA_merged
rm(list = setdiff(ls(), c("dataAA", "dataICPMS")))


# Defining lists for future loops
# excluding method blank and quality control from list of sites
sampleSites <- unique(filter(dataICPMS, site != "MB", site != "")$site)
metalsAnalyzed <- unique(dataICPMS$metal)

# Initiate loop, filter out cal data
calICPMS <- NULL

for(uniqueMetal in metalsAnalyzed){
  cal <- dataICPMS %>%
    filter(type == "Cal1" | type == "Cal2" | type =="Cal3") %>%
    filter(metal == uniqueMetal) %>%
    select(concentration, cps, rsd)
  
  #Step 3.2, perform weighted linear reg, and pull out relevant info from model
  w <- 1/ (cal$cps * cal$rsd)^2
  model <- lm(cal$cps ~ cal$concentration, weights = w)
  
  slope <- model$coefficients[2]
  intercept <- model$coefficients[1]
  slopeStd <- summary(model)$coefficients[2,2]
  interceptStd <- summary(model)$coefficients[1,2]
  
  plot(cal$cps ~ cal$concentration,
       xlab = paste("Concentration of ", uniqueMetal, "(ppb)"),
       ylab = "Counts per second")+
    abline(model, col = "red")+
    title(paste("Calibration for", uniqueMetal))
  
  equation <- tibble(metal = uniqueMetal, slope, slopeStd, intercept, interceptStd)
  calICPMS <- rbind(calICPMS, equation)
}

remove(equation, cal, slope, slopeStd, intercept, interceptStd, w, model, uniqueMetal)


