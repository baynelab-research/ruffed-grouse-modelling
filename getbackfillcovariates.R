#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(wildRtrax) #to download data from wildtrax
library(lubridate)

#2. Set root path for data on google drive----

root <- "G:/Shared drives/RUGR_LAB_PROJECT/DATA"

#A. DOWNLOAD DATA FROM WILDTRAX#######################

#1. Login to WildTrax----
config <- "login.R"
source(config)

wt_auth()

#2. Download summary report----
raw <- wt_download_report(project_id = 1321, sensor_id = "ARU", weather_cols = F, report = "summary")

write.csv(raw, file.path(root, "RUGR - BU LAB PROJECT - wildTrax summary report.csv"), row.names=FALSE)

#B. WRANGLE####

#1. Summarize RUGR data to visit covs----
visit <- raw %>%
    dplyr::select(location, recording_date, longitude, latitude, method, status, observer) %>%
    unique() %>%
    rename(date = recording_date) %>%
    mutate(date = ymd_hms(date),
           year = year(date))

#2. Summarize bird data per visit----
detection <- raw %>%
    rename(date = recording_date, species = species_code) %>%
    dplyr::filter(species=="RUGR") %>%
    mutate(date = ymd_hms(date)) %>%
    group_by(location, date) %>%
    summarize(abundance = sum(abundance)) %>%
    ungroup() %>%
    right_join(visit %>%
                   dplyr::select(location, date)) %>%
    mutate(abundance = ifelse(is.na(abundance), 0, abundance))


#C. ADD COVARIATE DATA####


#1. Load----
load("G:/My Drive/ABMI/Projects/BirdModels/Data/2wrangled.Rdata")


#3. Put covariates together----
covs <- dd %>%
    rename(location=SS, year = YEAR, equipment = EQUIP) %>%
    dplyr::select(location, year, AHM, FFP, MAP, MAT, MCMT, MWMT, PET, pAspen, NRNAME, NSRNAME, LUF_NAME, equipment) %>%
    right_join(visit)


#D. SAVE####
save(detection, covariate, file.path(root, "RUGR - BU LAB PROJECT - detections & ABMI covariates.Rdata"))
