#----------------- NUTRIENT EFFECT IN S:H RATIO and YII IN SSIDEREA

getwd()
source("STEPoneFunction.R")# R.Cunning steponeR function
library(plyr)
library(reshape2)
library(ggplot2)


      ##----------------------------------------------------------##
                           # 1. qPCR RATIOS 
      ##----------------------------------------------------------##
                          

# # # Get the raw data for Siderastrea siderea applying R.Cunning steponeR function----------


# Get list of plate files to read in
  Ssid.plates <- list.files(path="data/qPCR", pattern="csv$", full.names=T)
  Ssid.plates


# # # Run stepone function to get Ratios
  Ssid.Out <- steponeR(files=Ssid.plates, target.ratios=c("A.Ssid", "B.Ssid", "C.Ssid", "D.Ssid"), 
                     fluor.norm=list(A=-0.064, B=4.197, C=3.798, D=0, Ssid=7.416),
                     copy.number=list(A=1, B=1, C=50, D=3, Ssid=1),
                     ploidy=list(A=1, B=1, C=1, D=1, Ssid=2),
                     extract=list(A=0.813, B=0.813, C=0.813, D=0.813, Ssid=0.982))

# Target ratio results
  Ssid<-Ssid.Out$result

##----------------------------------------------------------##
# Labaling and Factors

# Parse sample names, times, treatments, colonies, plates
    sample.treatments <- rbind.fill(lapply(strsplit(as.character(Ssid$Sample.Name), split="_"), 
                                           function(X) data.frame(t(X))))
    colnames(sample.treatments) <- c("Treatment", "Time","Repli", "Spp","Sample")
    Ssid <- cbind(sample.treatments, Ssid[,-1])
    
    Sample <- rbind.fill(lapply(strsplit(as.character(Ssid$Sample), split="-"), 
                                function(X) data.frame(t(X))))
    colnames(Sample) <- c("Colony", "Core")
    Ssid <- cbind(Sample, Ssid)

# Keep unique sample ID to check reruns  
    Ssid$Sample.Time<-paste(Ssid$Time,Ssid$Sample, sep="_")
    Ssid$ID<-paste(Ssid$Sample.Time,Ssid$File.Name, sep=".")

##----------------------------------------------------------##

# Calculate total S/H ratio and D/C ratio and propD

# 0. If Clade only detected in one technical replicate, set its ratio to NA (becomes zero)
    Ssid$A.Ssid[which(Ssid$A.reps==1)] <- NA
    Ssid$B.Ssid[which(Ssid$B.reps==1)] <- NA
    Ssid$C.Ssid[which(Ssid$C.reps==1)] <- NA
    Ssid$D.Ssid[which(Ssid$D.reps==1)] <- NA

# 1. Rename cols and make NA=0
    colnames(Ssid)[which(colnames(Ssid) %in% c("A.Ssid", "B.Ssid","C.Ssid", "D.Ssid" ))] <- c("A.SH", "B.SH", "C.SH", "D.SH")  
    
    Ssid$A.SH[is.na(Ssid$A.SH)] <- 0
    Ssid$B.SH[is.na(Ssid$B.SH)] <- 0
    Ssid$C.SH[is.na(Ssid$C.SH)] <- 0
    Ssid$D.SH[is.na(Ssid$D.SH)] <- 0

# Get the ratios and log 10
  # Total ratio
      Ssid$tot.SH <- Ssid$A.SH + Ssid$B.SH + Ssid$C.SH + Ssid$D.SH  # Add A, B C and D to get total SH
      Ssid$logTot.SH <- log10(Ssid$tot.SH )  # Calculate log10 SH ratio
  # B ratio
      Ssid$logB.SH <- log10(Ssid$B.SH)  
  # C ratio
      Ssid$logC.SH <- log10(Ssid$C.SH)
  # D ratio
      Ssid$logD.SH <- log10(Ssid$D.SH)  

# Clade Proportion
  # D Proportion
    Ssid$D.Prp<-(Ssid$D.SH/Ssid$tot.SH)
  # C Proportion
    Ssid$C.Prp<-(Ssid$C.SH/Ssid$tot.SH)
  # B Proportion
    Ssid$B.Prp<-(Ssid$B.SH/Ssid$tot.SH)

    Ssid$logTot.SH[which(Ssid$tot.SH==0)] <- NA
    Ssid$logA.SH[which(Ssid$A.SH==0)] <- NA
    Ssid$logB.SH[which(Ssid$B.SH==0)] <- NA
    Ssid$logC.SH[which(Ssid$C.SH==0)] <- NA
    Ssid$logD.SH[which(Ssid$D.SH==0)] <- NA
 
 # Core clasification by clade Proportion 
    Ssid$Community[which(Ssid$D.Prp>=0.9)] <- "D"
    Ssid$Community[which(Ssid$D.Prp<=0.1)] <- "C"
    Ssid$Community[which(Ssid$D.Prp<0.9 & Ssid$D.Prp>0.1)] <- "DC"

  Ssid$Community<-factor(as.character(Ssid$Community), levels=c("C","DC","D"))
  
  propD <- Ssid$D.Prp
  hist(propD)
  # Almost all samples D dominated (>90% D)
 
  propDC <- Ssid$D.Prp[which(Ssid$D.Prp > 0 & Ssid$D.Prp < 0.9)]
  hist(propDC)
  range(propDC)


    ##----------------------------------------------------------##
                     # 2. qPCR DATA CLEANING 
    ##----------------------------------------------------------##

# 1. Check and remove NTC wells
    ntc <- Ssid[which(Ssid$Sample=="S0-NTC"), ]
    if(any(!is.na(ntc$CT))) warning("Template detected in NTC: interpret data with caution")
    Ssid <- droplevels(Ssid[!rownames(Ssid) %in% rownames(ntc), ])

# 2. Chose bw samples ran more than once
    ReRunA <- Ssid[duplicated(Ssid$Sample.Time),]
    
    n_RunA <- data.frame(table(Ssid$Sample.Time))
    colnames(n_RunA)<-c("Sample.Time","RanA")
    Ssid<-join(Ssid, n_RunA, type = "left")
    
    DuplicatesA <- Ssid[(Ssid$RanA>1),]
    write.csv(DuplicatesA, file = 'DuplicatesA.csv')

# 3. Remove bad replicates
    ToRem1<-read.csv("data/BadReplicates2.csv") # 11/25/15
    Ssid<-Ssid[!(Ssid$ID %in% ToRem1$ID),]
    
    n_RunB <- data.frame(table(Ssid$Sample.Time))
    ReRunB <- Ssid[duplicated(Ssid$Sample.Time),]
    colnames(n_RunB)<-c("Sample.Time","RanB")
    Ssid<-join(Ssid, n_RunB, type = "left")

# 4. Check and remove NoHSratio samples
    NoHSratio <- Ssid[which(Ssid$tot.SH==0), ]
    Ssid <- droplevels(Ssid[!rownames(Ssid) %in% rownames(NoHSratio), ])

    n_RunC <- data.frame(table(Ssid$Sample.Time))
    ReRunC <- Ssid[duplicated(Ssid$Sample.Time),]
    colnames(n_RunC)<-c("Sample.Time","RanC")
    Ssid<-join(Ssid, n_RunC, type = "left")

# List of dupplicated samples, should have 0 rows now
  DuplicatesC <- Ssid[(Ssid$RanC>1),]
  # write.csv(DuplicatesC, file = 'DuplicatesC.csv')

# 5. Check and remove samples with high SD
    # StDe2S <- Ssid[which(Ssid$Ssid.CT.sd>2), ]
    # StDe1.5S <- Ssid[which(Ssid$Ssid.CT.sd>1.5), ]
    StDe2all <- Ssid[which((Ssid$C.CT.sd)>1.5|(Ssid$D.CT.sd>2)|(Ssid$Ssid.CT.sd>2)), ]
    # StDe1.5all <- Ssid[which((Ssid$C.CT.sd)>1.5 |(Ssid$D.CT.sd>1.5) |(Ssid$Ssid.CT.sd>1.5)), ]
    #Ssid <- droplevels(Ssid[!rownames(Ssid) %in% rownames(StDe2all), ])
    # write.csv(StDe2all, file = 'StDe2all.csv')

summary(Ssid)

# Export data if want local backup
  # write.csv(Ssid, file = 'Ssid.csv')

# End of qPCR data cleaning



        ##----------------------------------------------------------##
            # 3. Explore SH and Clade proportions with ggplot
        ##----------------------------------------------------------##
  
#library(ggplot2)

# Re order factors
  
  Ssid$Treatment <- as.character(Ssid$Treatment)
  Ssid$Treatment[Ssid$Treatment == "C"] <- "Control"
  Ssid$Treatment[Ssid$Treatment == "D"] <- "Dark"
  
  Ssid$Treatment <- factor (as.character(Ssid$Treatment), levels=c("CO2","Control","Fe", "N", "NP", "NPF","Dark"))
    # print(levels(Ssid$Treatment))
    
  Ssid$Time <- factor (as.character(Ssid$Time), levels=c("T0","T1","T2","T3"))

  Ssid$Colony<-factor(as.character(Ssid$Colony), 
                    levels=c("S1","S4","S3","S5","S6","S2"))  # D dominance order

# Graphs with GGPLOT
    
  # D proportion frecuency / time or treatment
    #Community <- ggplot(data=Ssid, aes(Ssid$D.Prp)) 
    #Community + geom_histogram(binwidth = 0.1) + facet_grid(Time~., margin=T)
    #Community + geom_histogram(binwidth = 0.1) + facet_grid(Treatment~., margin=T)

    #logSHRTime <- ggplot(Ssid, aes(factor(Time), logTot.SH))
    #logSHRTime + geom_boxplot(aes(fill=factor(Treatment)))+ facet_grid(Colony~., margins=TRUE)
  
  # Total SH by Treatmemnt*Time*Colony  
    logSHRTreatment <- ggplot(Ssid, aes(factor(Treatment), logTot.SH))
    # logSHRTreatment + geom_boxplot(aes(fill=factor(Time)))
    logSHRTreatment + geom_boxplot(aes(fill=factor(Time)))+ facet_grid(Colony~., margins=TRUE)
    
    #     logSHRcomm<- ggplot(Ssid, aes(factor(Time), logTot.SH))
    #     logSHRcomm + geom_boxplot(aes(fill=factor(Treatment)))
    #     logSHRcomm + geom_boxplot(aes(fill=factor(Community)))
        
    #     logSHRcomm2<- ggplot(Ssid, aes(factor(Treatment), logTot.SH))
    #     logSHRcomm2 + geom_boxplot(aes(fill=factor(Community)))+ facet_grid(Time~., margins=TRUE)
    
  # Change in log C by Treatmemnt*Time*Colony
     SHRC <- ggplot(Ssid, aes(factor(Treatment), logC.SH))
     # SHRC + geom_boxplot(aes(fill=factor(Time)))
     SHRC + geom_boxplot(aes(fill=factor(Time)))+ facet_grid(Colony~., margins=TRUE)
    # CH increased in Dark, N, Fe?, but not NP, or NPF
  
# Change in log D by Treatmemnt*Time*Colony
     SHRD <- ggplot(Ssid, aes(factor(Treatment), logD.SH))
     SHRD + geom_boxplot(aes(fill=factor(Time)))+ facet_grid(Colony~., margins=TRUE)
    # Did D increased in Dark?

  # Change in D proportion by Treatmemnt*Time*Colony 
     Dprop <- ggplot(Ssid, aes(factor(Treatment), D.Prp))
     Dprop + geom_boxplot(aes(fill=factor(Time)))+ facet_grid(Colony~., margin=T)

  # Change in C proportion by Treatmemnt*Time*Colony  
    Cprop <- ggplot(Ssid, aes(factor(Treatment), C.Prp))
    Cprop + geom_boxplot(aes(fill=factor(Time)))+ facet_grid(Colony~., margin=T) 
  # Colony S2 does not change at all!!! Completely D dominated


      ##----------------------------------------------------------##
                # 4. Join IPAM data and qPCR DATA
      ##----------------------------------------------------------##

qPCR <- Ssid[,c("Sample.Time", "Colony","Community", "Treatment","Time","Repli","Spp","Core","C.SH","D.SH","tot.SH", "logC.SH","logD.SH", "logTot.SH","D.Prp")]
IPAM_Final<- read.csv("data/IPAM-data.csv")

All.data<-join(qPCR, IPAM_Final, by = "Sample.Time", type = "inner", match = "all")

# Remove duplicated variables
  All.data$Core <- NULL
  All.data$Spp_ <- NULL
  All.data$Treatment_ <- NULL
  All.data$Rep_ <- NULL

# Keep a local file
# write.csv(All.data, file = 'AllData.csv')


##----------------------------------------------------------##
                      # 5. DATA ANALISIS
##----------------------------------------------------------##

# LMER Models for different set of data

# Libraries
  library (lattice)
  library(lme4)
  library(lmerTest)
  library(effects)

# Histograms
  NorTot <- qplot(logTot.SH, data=All.data, geom="histogram")
  NorD <-  qplot(logD.SH, data=All.data, geom="histogram")
  NorC <- qplot(logC.SH, data=All.data, geom="histogram")
  NorYII <-qplot(Y.II., data=All.data, geom="histogram")
  
  # NorDProp <-qplot(D.Prp, data=All.data, geom="histogram")
  # All.data$D.Prp.t <- asin(sqrt(All.data$D.Prp))
  # NorDProp.t <-qplot(D.Prp.t, data=All.data, geom="histogram")

# -------------------------------------------------------------
                      # A. ALL TREATMENTS 

#Tot SH
  Mod_Tot_All <- lmer(logTot.SH ~ Treatment * Time + (Time|Colony), data=All.data)
  step(Mod_Tot_All)
  anova(Mod_Tot_All)
  plot(Effect(c("Time", "Treatment"), Mod_Tot_All), x.var="Time", multiline=T, ci.style="bars")
  summary(Mod_Tot_All)

#C:H
  Mod_C_All <- lmer(logC.SH ~ Treatment * Time + (Time|Colony), data=All.data)
  step(Mod_C_All)
  anova(Mod_C_All)
  plot(Effect(c("Time", "Treatment"), Mod_C_All), x.var="Time", multiline=T, ci.style="bars")
  summary(Mod_C_All)

#D:H
  Mod_D_All <- lmer(logD.SH ~ Treatment * Time + (Time|Colony), data=All.data)
  step(Mod_D_All)
  anova(Mod_D_All)
  plot(Effect(c("Time", "Treatment"), Mod_D_All), x.var="Time", multiline=T, ci.style="bars")
  summary(Mod_D_All)

#YII
  Mod_YII_All <- lmer(Y.II. ~ Treatment * Time + (Time|Colony), data=All.data)
  step(Mod_YII_All)
  anova(Mod_YII_All)
  plot(Effect(c("Time", "Treatment"), Mod_YII_All), x.var="Time", multiline=T, ci.style="bars")
  summary(Mod_YII_All)

# -------------------------------------------------------------
                  # B. NO CO2 AND DARK

