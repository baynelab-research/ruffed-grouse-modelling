#NOTES


#PREAMBLE####

#1. Load packages----
library(tidyverse)
library(QPAD)
library(maptools)
library(intrval)
library(raster)
library(lme4)
library(MASS)

#2. Set google drive root----
root <- "G:/Shared drives/RUGR_LAB_PROJECT"

#3. Load data----
load(file.path(root, "DATA", "RUGR - BU LAB PROJECT - detections & ABMI covariates.Rdata"))

#4. Read in raw data----
raw <- read.csv(file.path(root, "DATA", "RUGR - BU LAB PROJECT - wildTrax summary report - 2023-02-08.csv"))

#GET OFFSETS####

#1. Load offsets----
load_BAM_QPAD(3)

#2. Join with covariates & format----
dat.x <- covariate %>%
    dplyr::select(location, year, date, latitude, longitude, pAspen, MAP) %>%
    left_join(detection) %>%
    mutate(datetime = ymd_hms(date),
           date = str_sub(as.character(datetime), 1, 10),
           time = str_sub(as.character(datetime), -8, -1),
           dur = 1,
           dis = Inf,
           tagmeth = "ARU") %>%
    dplyr::rename(lat = latitude, lon = longitude) %>%
    dplyr::filter(!is.na(lat),
                  !is.na(pAspen))

#3. Set up qpad offsets stuff----
setwd("C:/Users/elly/Documents/BAM/QPAD/qpad-offsets")

source("functions.R")

#4. Read in rasters----
rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)

#5. Make x----
x <- make_x(dat.x, tz="local")

#6. Calculate offsets-----
off <- make_off("RUGR", x, useMeth = "n")

#7. Put together----
dat <- cbind(dat.x, off, x)

#8. Reset wd----
setwd("C:/Users/elly/Documents/BayneLab/ruffed-grouse-modelling/qpad")

#SINGLE VISIT QPAD####

#1. Set seed----
set.seed(1234)

#2. Randomly select one visit per location----
dat1 <- dat %>%
    group_by(location) %>%
    sample_n(1) %>%
    ungroup()

#3. Visualize----
ggplot(dat1) +
#    geom_point(aes(x=pAspen, y=abundance)) +
    geom_smooth(aes(x=pAspen, y=abundance))

#4. Model----
m1 <- glm(abundance ~ pAspen, offset=offset, data = dat1, family="poisson")
summary(m1)

#5. Predict----
newdat <- data.frame(pAspen = seq(0, 1, 0.01))

pred1 <- data.frame(predict(m1, newdat, se=TRUE, type="response")) %>%
    cbind(newdat)

#6. Plot----
ggplot(pred1) +
    geom_ribbon(aes(x=pAspen, ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit), alpha = 0.5) +
    geom_line(aes(x=pAspen, y=fit))

#MULTIVISIT QPAD####

#1. Set seed----
set.seed(1234)

#2. Visualize----
ggplot(dat) +
    #    geom_point(aes(x=pAspen, y=abundance)) +
    geom_smooth(aes(x=pAspen, y=abundance))

#4. Model----
mall <- glmer(abundance ~ pAspen + (1|location), offset=offset, data = dat, family="poisson")
summary(mall)

#5. Predict----
newdat <- data.frame(pAspen = seq(0, 1, 0.01))

predall <- data.frame(fit = predict(mall, newdat, type="response", re.form = NA)) %>%
    cbind(newdat)

ciall <- confint(mall)

#6. Plot----
ggplot(pred1) +
    geom_ribbon(aes(x=pAspen, ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit), alpha = 0.5) +
    geom_line(aes(x=pAspen, y=fit)) +
    geom_line(data=predall, aes(x=pAspen, y=fit), colour="blue")

#MAXVISIT QPAD####

#1. Set seed----
set.seed(1234)

#2. Randomly select one visit per location that has max abundance----
datmax <- dat %>%
    group_by(location) %>%
    mutate(total = max(abundance)) %>%
    dplyr::filter(abundance == total) %>%
    sample_n(1) %>%
    ungroup()

#3. Visualize----
ggplot(datmax) +
    #    geom_point(aes(x=pAspen, y=abundance)) +
    geom_smooth(aes(x=pAspen, y=abundance))

#4. Model----
mmax <- glm(abundance ~ pAspen, offset=offset, data = datmax, family="poisson")
summary(mmax)

#5. Predict----
newdat <- data.frame(pAspen = seq(0, 1, 0.01))

predmax <- data.frame(predict(mmax, newdat, se=TRUE, type="response")) %>%
    cbind(newdat)

#6. Plot----
ggplot(pred1) +
    geom_ribbon(aes(x=pAspen, ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit), alpha = 0.3) +
    geom_ribbon(data=predmax, aes(x=pAspen, ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit), alpha = 0.3, fill="red") +
    geom_line(aes(x=pAspen, y=fit)) +
    geom_line(data=predall, aes(x=pAspen, y=fit), colour="blue") +
    geom_line(data=predmax, aes(x=pAspen, y=fit), colour="red")

ggsave(filename = file.path(root, "FIGURES", "QPADPredictions-paspen.jpeg"), width = 6, height = 4)

#COMPARE###
ests <- data.frame(summary(m1)[["coefficients"]]) %>%
    mutate(model = "Singlevisit",
           variable = c("Intercept", "pAspen")) %>%
    rbind(data.frame(summary(mmax)[["coefficients"]]) %>%
              mutate(model = "Maxvisit",
                     variable = c("Intercept", "pAspen"))) %>%
    rbind(data.frame(summary(mall)[["coefficients"]]) %>%
              mutate(model = "Multivisit",
                     variable = c("Intercept", "pAspen")))

ggplot(ests) +
    geom_point(aes(x=variable, y=Estimate, colour=model), position = position_dodge(width = 1)) +
    geom_errorbar(aes(x=variable, ymin=Estimate-1.96*Std..Error, ymax=Estimate+1.96*Std..Error, colour=model), position = position_dodge(width = 1))

ggsave(filename = file.path(root, "FIGURES", "QPADCoefficients-paspen.jpeg"), width = 6, height = 4)



