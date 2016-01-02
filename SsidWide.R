getwd()
library(plyr)
library(MASS)
library(effects)
library(reshape2)
library(ggplot2)


# -------------------------------------------------------------------------------------------------
# ANALYSES

# Modify Ssid.csv and save it as data.csv

# -------------------------------------------------------------------------------------------------

getwd()
data<-read.csv("data.csv")
library(tidyr)
qPCRwide<-spread(data, Time, measurement)
data_wide

data$ID <- as.factor(data$ID)
data$Colony <- as.factor(data$Colony)
data$Treatment <- as.factor(data$Treatment)
data$Time <- as.factor(data$Time)
data$Repli <- as.factor(data$Repli)
data$Community<- as.factor(data$Community)

data$Treatment <- factor (as.character(data$Treatment), levels=c("CO2","C","Fe", "N", "NP", "NPF","D"))
print(levels(data$Treatment))

