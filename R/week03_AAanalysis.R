## Week 03 Data Analysis for AA data
## Gillian McGinnis
## Created 06 November 2020
## Updated 11 November 2020

library(tidyverse)
library(ggthemes)
source("R/week02_fullTidy.R")
rm(list = setdiff(ls(), c("AA_merged", "dataKey_tidy")))

## Adjusting dataframes
data_aa <- AA_merged %>%
  mutate(xRsd = as.numeric(as.character(case_when(xRsd != "HIGH" ~ xRsd)))) %>%
  clean_names() %>%
  mutate(rsd = x_rsd)

sample_key <- dataKey_tidy %>%
  clean_names()

## Cleaning environment
rm(list = setdiff(ls(), c("data_aa", "sample_key")))

## Creating a list of sample sites by excluding NAs and Method Blanks
sample_sites <- unique(filter(data_aa, site != "MB", site != "")$site)

## Code based on ICPMS analysis
# Initiate loop, filter out cal data
cal_aa <- NULL
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
  
equation <- tibble(slope, slopeStd, intercept, interceptStd)
#cal_aa <- rbind(cal_aa, equation)
cal_aa <- rbind(equation)
cal_aa


## Alternatively, here's my way.
myResults <- summary(lm(mean_abs ~ concentration, weights = w, data = cal))

AA_cal_curve <- ggplot(cal, aes(x = concentration, y = mean_abs))+
  geom_smooth(method = "lm", se = FALSE, color = "red")+
  geom_point(shape = 1)+
  theme_tufte()+
  theme(axis.line = element_line("black"),
        text = element_text("Times", size=12))+
  labs(x = "Concentration of Cr (ppb)", y = "Counts per second", title = "Calibration curve for Cr")

ggsave("AA_cal_curve.png", plot = AA_cal_curve, path = "figures/")

## Environment clean
remove(equation, cal, slope, slopeStd, intercept, interceptStd, w, model)

cal_aa <- cal_aa %>%
  clean_names()

## Step 4.1: creating a function to analyze the samples
#inputs: uniqueSite (as a character)
#outputs: concentration vector
sample_analysis <- function(unique_site){
#unique_site <- "A"  
concentration_data <- NULL
#  for(unique_metal in metals_analyzed){
  sample <- filter(data_aa,
                     #metal == unique_metal,
                   site == unique_site)
  data <- NULL
  for(ID in sample$sample_key){
    sample_data <- filter(sample, sample_key == ID)
    #cal <- filter(calAA, metal == unique_metal)
    cal <- cal_aa
      
    m <- cal$slope
    b <- cal$intercept
    y <- sample_data$mean_abs
      
    b_e <- cal$intercept_std
    m_e <- cal$slope_std
      
    x <- (y-b)/m
      
      #RSD <- sample_data$rsd
    RSD <- ((sample_data$rsd/100)*sample_data$mean_abs)
    #CPS <- sample_data$cps
    ABS <- sample_data$mean_abs
      
    e_yb <- sqrt((RSD)^2 + (b_e)^2)
    #yb <- CPS-b
    yb <- ABS-b
    e_x <- x*sqrt((e_yb/yb)^2 + (m_e/m)^2)
      
    data <- rbind(data, data_frame(sample_key = ID, x, e_x))
      
    if(unique_site != "MB"){
      concentration_data <- data_frame(sample_key = sample_data$sample_key,
                                       analyst = sample_data$analyst,
                                       #metal = unique_metal,
                                       site = unique_site,
                                       conc_dil = x,
                                       conc_dil_error = e_x) %>%
        rbind(concentration_data)
      }
    }
  if(unique_site == "MB"){
      x <- mean(data$x)
      e_x <- sd(data$x)
      concentration_data <- data_frame(site = unique_site,
                                       #metal = unique_metal,
                                       conc_dil = x,
                                       conc_dil_error = e_x) %>%
        rbind(concentration_data)
    }
  #} #closing bracket for unique metal if statement
  return(concentration_data)
}

## Step 5: create a function that runs a diff function on ea of the soil sample sites
# inputs: a function
# outputs: a dataframe with the function outputs from each site
run_sites <- function(Function){
  value <- NULL
  for(sites in sample_sites){
    site_value <- Function(sites)
    value <- rbind(site_value, value)
  }
  return(value)
}

MB <- sample_analysis("MB")
uncor_sample <- run_sites(sample_analysis)

MB
uncor_sample


## Step 7: correct for the MB and perform error prop as needed
sample_data_mb <- NULL

conc_dil_blanked <- uncor_sample$conc_dil-MB$conc_dil
conc_dil_blanked_error <- sqrt(uncor_sample$conc_dil_error)^2 + (MB$conc_dil_error)^2

sample_data_mb <- uncor_sample %>%
  mutate(conc_dil_blanked, conc_dil_blanked_error) %>%
  rbind(sample_data_mb)

sample_data_mb


## step 8: define dilution factors and measurement errors
vol_e <- 1
mass_e <- 0.001
dil_1010_e <- sqrt(1^2 + 10^2)
dil_e <- sqrt((dil_1010_e/1010)^2 + (1/10)^2) ## error in 101 dilution factor

sample_data <- merge(data_aa, sample_data_mb) %>%
  unique()%>%
  mutate(conc_blanked = conc_dil_blanked * (total_volume/1000)/(mass_of_soil/1000),
         conc_blanked_error = conc_blanked*
           sqrt((conc_dil_blanked_error/conc_dil_blanked)^2+
                  (dil_e/101)^2 +
                  (mass_e/mass_of_soil)^2 +
                  (vol_e/total_volume)^2),
         conc_unblanked = conc_dil*(total_volume/1000)/(mass_of_soil/1000),
         conc_unblanked_error = conc_unblanked*
           sqrt((conc_dil_error/conc_dil)^2+
                  (dil_e/101)^2 +
                  (mass_e/mass_of_soil)^2 +
                  (vol_e/total_volume)^2)) %>%
  select(!c(concentration, type, mass_of_soil, total_volume, mean_abs, rsd, conc_dil_blanked, conc_dil_blanked_error, conc_dil, conc_dil_error))

# purging the environment
rm(list = ls()[!ls() %in% c("data_aa", "sample_data")])