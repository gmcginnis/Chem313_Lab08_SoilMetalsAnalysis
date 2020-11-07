## Week 03 Data Analysis for ICPMS data
## Gillian McGinnis
## Created 06 November 2020
## Updated 06 November 2020

source("R/week02_fullTidy.R")
rm(list = setdiff(ls(), c("AA_merged", "ICPMS_merged", "dataKey_tidy")))

dataICPMS <- ICPMS_merged %>%
  clean_names("small_camel")
dataAA <- AA_merged
dataKey <- dataKey_tidy
rm(list = setdiff(ls(), c("dataAA", "dataICPMS", "dataKey")))


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


## I prefer camelCase but that seems to be causing errors. Over to snake_case we go. RIP my environment, you had a good run.
ICPMS <- dataICPMS %>%
  clean_names()
sample_key <- dataKey %>%
  clean_names()
ICPMS_cal <- calICPMS %>%
  clean_names()
metals_analyzed <- metalsAnalyzed
sample_sites <- sampleSites
# Step 4.1: creating a function to analyze the samples
#inputs: uniqueSite (as a character)
#outputs: concentration vector
sample_analysis <- function(unique_site){
  concentration_data <- NULL
  for(unique_metal in metals_analyzed){
    sample <- filter(ICPMS, metal == unique_metal, site == unique_site)
    data <- NULL
    
    for(ID in sample$sample_key){
      sample_data <- filter(sample, sample_key == ID)
      cal <- filter(ICPMS_cal, metal == unique_metal)
      
      m <- cal$slope
      b <- cal$intercept
      y <- sample_data$cps
      
      b_e <- cal$intercept_std
      m_e <- cal$slope_std
      
      x <- (y-b)/m
      
      RSD <- sample_data$rsd
      CPS <- sample_data$cps
      
      e_yb <- sqrt((RSD)^2 + (b_e)^2)
      yb <- CPS-b
      e_x <- x*sqrt((e_yb/yb)^2 + (m_e/m)^2)
      
      data <- rbind(data, data_frame(sample_key = ID, x, e_x))
      
      if(unique_site != "MB"){
        concentration_data <- data_frame(sample_key = sample_data$sample_key,
                                         analyst = sample_data$analyst,
                                         metal = unique_metal,
                                         site = unique_site,
                                         conc_dil = x,
                                         conc_dil_error = e_x) %>%
          rbind(concentration_data)
      }
    }
    if(unique_site == "MB"){
      x <- mean(data$x)
      e_x <- sd(data$x)
      concentration_data <- data_frame(metal = unique_metal,
                                       site = unique_site,
                                       conc_dil = x,
                                       conc_dil_error = e_x) %>%
        rbind(concentration_data)
    }
  }
  return(concentration_data)
}

MB <- sample_analysis("MB")
uncor_sample <- runSites(sample_analysis)

MB
uncor_sample


sample_data_mb <- NULL

for(unique_metal in metals_analyzed){
  MB_metal <- filter(MB, metal == unique_metal)
  sample_metal <- filter(uncor_sample, metal == unique_metal)
  conc_dil_blanked <- sample_metal$conc_dil-MB_metal$conc_dil
  
  conc_dil_blanked_error <- sqrt(sample_metal$conc_dil_error)^2 + (MB_metal$conc_dil_error)^2
  
  sample_data_mb <- sample_metal %>%
    mutate(conc_dil_blanked, conc_dil_blanked_error) %>%
    rbind(sample_data_mb)
}

sample_data_mb


#step 8
vol_e <- 1
mass_e <- 0.001
dil_1010_e <- sqrt(1^2 + 10^2)
dil_e <- sqrt((dil_1010_e/1010)^2 + (1/10)^2) ## error in 101 dilution factor

sample_data <- merge(ICPMS, sample_data_mb) %>%
  unique()%>%
  mutate(conc_blanked = conc_dil_blanked * (total_volume/1000)/(mass_of_soil/1000)*101,
         conc_blanked_error = conc_blanked*
           sqrt((conc_dil_blanked_error/conc_dil_blanked)^2+
                  (dil_e/101)^2 +
                  (mass_e/mass_of_soil)^2 +
                  (vol_e/total_volume)^2)) %>%
  select(!c(concentration, type, mass_of_soil, total_volume, cps, rsd, conc_dil_blanked, conc_dil_blanked_error, conc_dil, conc_dil_error))

# purging the environment
rm(list = ls()[!ls() %in% c("ICPMS", "sample_data")])
