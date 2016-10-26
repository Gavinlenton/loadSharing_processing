# Basic data analysis for ROM data

# set working directory
setwd("~/Google Drive/LS_main_data_collection/ElaboratedData/ROM/")
hip_flexion_ROM <- read.csv("hip_flexion_ROM.csv", header = TRUE, sep = ",", na.strings = "NaN")
shoulder_FF_ROM <- read.csv("shoulder_FF_ROM.csv", header = TRUE, sep = ",", na.strings = "NaN")
UUA_ROM <- read.csv("UUA_ROM.csv", header = TRUE, sep = ",", na.strings = "NaN")
trunk_flexion_ROM <- read.csv("trunk_flexion_ROM.csv", header = TRUE, sep = ",", na.strings = "NaN")

# Check all variables for outliers
summary(hip_flexion_ROM)
summary(trunk_flexion_ROM)
summary(shoulder_FF_ROM)
# UUA looks good
summary(UUA_ROM)

# Convert super low values to NA
hip_flexion_ROM[,2:14][hip_flexion_ROM[,2:14] < 55] <- NA
shoulder_FF_ROM[,2:14][shoulder_FF_ROM[,2:14] < 90] <- NA
UUA_ROM[,2:14][UUA_ROM[,2:14] < 90] <- NA
trunk_flexion_ROM[,2:14][trunk_flexion_ROM[,2:14] < 70] <- NA

# Create subset of the data for 15 kg
hip_flexion_ROM.sub15kg <- subset(hip_flexion_ROM, select=c("TBAS15","CRYE15","TYR15","USMC15","CORE15","SORD15"))
shoulder_FF_ROM.sub15kg <- subset(shoulder_FF_ROM, select=c("TBAS15","CRYE15","TYR15","USMC15","CORE15","SORD15"))
UUA_ROM.sub15kg <- subset(UUA_ROM, select=c("TBAS15","CRYE15","TYR15","USMC15","CORE15","SORD15"))
trunk_flexion_ROM.sub15kg <- subset(trunk_flexion_ROM, select=c("TBAS15","CRYE15","TYR15","USMC15","CORE15","SORD15"))

# Create subset of the data for 30 kg
hip_flexion_ROM.sub30kg <- subset(hip_flexion_ROM, select=c("TBAS30","CRYE30","TYR30","USMC30","CORE30","SORD30"))
shoulder_FF_ROM.sub30kg <- subset(shoulder_FF_ROM, select=c("TBAS30","CRYE30","TYR30","USMC30","CORE30","SORD30"))
UUA_ROM.sub30kg <- subset(UUA_ROM, select=c("TBAS30","CRYE30","TYR30","USMC30","CORE30","SORD30"))
trunk_flexion_ROM.sub30kg <- subset(trunk_flexion_ROM, select=c("TBAS30","CRYE30","TYR30","USMC30","CORE30","SORD30"))

# Summary of 15 kg data
summary_HF_15kg <- dfapply(hip_flexion_ROM.sub15kg, favstats)
summary_shoulderFF_15kg <- dfapply(shoulder_FF_ROM.sub15kg, favstats)
summary_UUA_15kg <- dfapply(UUA_ROM.sub15kg, favstats)
summary_TF_15kg <- dfapply(trunk_flexion_ROM.sub15kg, favstats)

# Summary of 30 kg data
summary_HF_30kg <- dfapply(hip_flexion_ROM.sub30kg, favstats)
summary_shoulderFF_30kg <- dfapply(shoulder_FF_ROM.sub30kg, favstats)
summary_UUA_30kg <- dfapply(UUA_ROM.sub30kg, favstats)
summary_TF_30kg <- dfapply(trunk_flexion_ROM.sub30kg, favstats)


