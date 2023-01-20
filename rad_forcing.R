#### RADIATIVE FORCING
rm(list=ls())
library(tidyr)
library(tibble)

## set working directory and read file
setwd("~/OneDrive - Western Washington University/Thesis/Manuscripts/Bagley UAV Paper/Radiative forcing")
rad = read.csv("radiative forcing.csv")

## subset data by date
rad_702 <- subset(rad, rad$date == "7/2/21")
rad_730 <- subset(rad, rad$date == "7/30/21")

## subset by wavelengths
rad_702 <- subset(rad_702, rad_702$wavelength > 340)
rad_702 <- subset(rad_702, rad_702$wavelength < 860)
rad_730 <- subset(rad_730, rad_730$wavelength > 340)
rad_730 <- subset(rad_730, rad_730$wavelength < 860)


## calculate instantaneus radiative forcing for each algae sample
irf.702 = data.frame(
  date = "2021-07-02",
  a1 = sum((rad_702$clean.reflectance-rad_702$algae.reflectance.1)*rad_702$global.irradiance*10),
  a2 = sum((rad_702$clean.reflectance-rad_702$algae.reflectance.2)*rad_702$global.irradiance*10),
  a3 = sum((rad_702$clean.reflectance-rad_702$algae.reflectance.3)*rad_702$global.irradiance*10),
  a4 = sum((rad_702$clean.reflectance-rad_702$algae.reflectance.4)*rad_702$global.irradiance*10),
  a5 = sum((rad_702$clean.reflectance-rad_702$algae.reflectance.5)*rad_702$global.irradiance*10),
  a6 = sum((rad_702$clean.reflectance-rad_702$algae.reflectance.6)*rad_702$global.irradiance*10),
  a7 = sum((rad_702$clean.reflectance-rad_702$algae.reflectance.7)*rad_702$global.irradiance*10),
  a8 = sum((rad_702$clean.reflectance-rad_702$algae.reflectance.8)*rad_702$global.irradiance*10),
  a9 = sum((rad_702$clean.reflectance-rad_702$algae.reflectance.9)*rad_702$global.irradiance*10)
)

irf.702.long = gather(irf.702, sample, irf, a1:a9, factor_key=TRUE) #transform data to long format

irf.730 = data.frame(
  date = "2021-07-30",
  a1 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.1)*rad_730$global.irradiance*10),
  a2 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.2)*rad_730$global.irradiance*10),
  a3 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.3)*rad_730$global.irradiance*10),
  a4 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.4)*rad_730$global.irradiance*10),
  a5 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.5)*rad_730$global.irradiance*10),
  a6 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.6)*rad_730$global.irradiance*10),
  a7 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.7)*rad_730$global.irradiance*10),
  a8 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.8)*rad_730$global.irradiance*10),
  a9 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.9)*rad_730$global.irradiance*10),
  a10 = sum((rad_730$clean.reflectance-rad_730$algae.reflectance.10)*rad_730$global.irradiance*10)
)

irf.730.long = gather(irf.730, sample, irf, a1:a10, factor_key=TRUE) #transform data to long format

## summary of IRF
summary(irf.702.long)
sd(irf.702.long$irf)
summary(irf.730.long)
sd(irf.730.long$irf)

## all data IRF
irf.all <- rbind(irf.702.long, irf.730.long)
summary(irf.all)

## 334,000 J needed to melt 1 kg of snow (Cohen 1994) -> XXX m3 of snow melted by algae
## 1 W = 1 J/s
## snow density = 600 kg/m3
## july 2, 2021 = 15:53 hr of daylight = 57180 s
## july 30, 2021 = 15:01 hr of daylight = 54060 s


## estimates of MF and m3 melt per day
irf.day.702 = data_frame(
  MJ.702 = (((irf.702.long$irf)*57180*1352.18)/1000000))
irf.day.702$melt.702 = ((irf.day.702$MJ.702*1000000)/(334000*600))

irf.day.730 = data_frame(
  MJ.730 = (((irf.730.long$irf)*54060*556.07)/1000000))
irf.day.730$melt.730 = ((irf.day.730$MJ.730*1000000)/(334000*600))

#get summary stats by group
summary(irf.day.702)
sd(irf.day.702$MJ.702)
sd(irf.day.702$melt.702)

summary(irf.day.730)
sd(irf.day.730$MJ.730)
sd(irf.day.730$melt.730)


## estimates per season (92 days)
## growing season from June 1 through August 31 -> 92 days

irf.out.702.season = irf.day.702*92
irf.out.730.season = irf.day.730*92

summary(irf.out.702.season)
summary(irf.out.730.season)


#error propogation
sd.702.1 <-(sd(irf.702.long$irf)*57180*1352.18)/1000000
sd.702.2 <- ((sd.702.1*1000000)/(334000*600))

sd.730.1 <-(sd(irf.730.long$irf)*54060*556.07)/1000000
sd.730.2 <- ((sd.730.1*1000000)/(334000*600))
