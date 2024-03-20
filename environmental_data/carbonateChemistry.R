library(lubridate) # for working with dates
library(ggplot2)  # for creating graphs
library(scales)   # to access breaks/formatting functions
library(gridExtra) # for arranging plots

library(readxl)
library(plyr)
library(dplyr)
library(ggpubr)
library(openair)
library(seacarb)

library(tidyverse)

############Data Import####################
#temperature
temperatureDataTable1 <- read_excel(paste("D:/googleDrive/remoteWork/papers/urbanCoral/temperature","/temperatureCol12.xlsx",sep=""))
temperatureDataTable2 <- read_excel(paste("D:/googleDrive/remoteWork/papers/urbanCoral/temperature","/temperatureCol3.xlsx",sep=""))
temperatureDataRaw <- rbind(temperatureDataTable1,temperatureDataTable2)
temperatureDataQAQC <- filter(temperatureDataRaw, good == "y" & fiveMin == "y")
temperatureDataQAQC$rDateTime <- ymd_hms(temperatureDataQAQC$dateTimeEST)
temperatureDataQAQC$rDateTime <- round_date(temperatureDataQAQC$rDateTime, unit = "minute")
temperatureDataQAQC$season <- factor(temperatureDataQAQC$season,levels = c("winter","spring","summer","autumn"))
tempEm <- filter(temperatureDataQAQC,  site == "em" )#site specific data used for coordinating with bottle data
tempNmac <- filter(temperatureDataQAQC,  site == "nmac" )#site specific data used for coordinating with bottle data
tempSmac <- filter(temperatureDataQAQC,  site == "smac" )#site specific data used for coordinating with bottle data
tempStar <- filter(temperatureDataQAQC,  site == "star" )#site specific data used for coordinating with bottle data
#seafet
phDataRaw <- read_excel(paste("D:/googleDrive/remoteWork/papers/urbanCoral/seafet","/seafetCombined.xlsx",sep=""))
phDataRaw$rDateTime <- ymd_hms(phDataRaw$dateTimeEST)#transform the date and time with lubridate
phDataRaw$rDateTime <- round_date(phDataRaw$rDateTime, unit = "minute")
numbRows = nrow(phDataRaw)
for(i in 1:numbRows){
  if (phDataRaw$rDateTime[i]>"2020-5-26" && phDataRaw$rDateTime[i]<"2020-10-15"){phDataRaw$precipitation[i] = "wet"}
  else if (phDataRaw$rDateTime[i]>"2019-5-26" && phDataRaw$rDateTime[i]<"2019-10-15"){phDataRaw$precipitation[i] = "wet"}
  else if (phDataRaw$rDateTime[i]>"2018-5-26" && phDataRaw$rDateTime[i]<"2018-10-15"){phDataRaw$precipitation[i] = "wet"}
  else {phDataRaw$precipitation[i] = "dry"}
}
phEm <- filter(phDataRaw,  site == "em" )#site specific data used for coordinating with bottle data
phNmac <- filter(phDataRaw,  site == "nmac" )#site specific data used for coordinating with bottle data
phSmac <- filter(phDataRaw,  site == "smac" )#site specific data used for coordinating with bottle data
phStar <- filter(phDataRaw,  site == "star" )#site specific data used for coordinating with bottle data
#bottle
bottleDataRaw <- read_excel(paste("D:/googleDrive/remoteWork/papers/urbanCoral/bottleSamples","/bottleICEQAQCWithAdditions.xlsx",sep=""))



################################################################################################################
#filter added for long bag holding time
bottleDataRaw <- filter(bottleDataRaw,  bagLessMonth == "y" )#remove samples that have sat in bags longer than a month


bottleDataRaw$rDateTime <- ymd_hms(bottleDataRaw$dateTimeEST)#transform the date and time with lubridate
bottleDataRaw$rDateTime30Rounded <- round_date(bottleDataRaw$rDateTime, "30 minutes")#transform the date and time with lubridate
numbRows = nrow(bottleDataRaw)
for(i in 1:numbRows){
  if (bottleDataRaw$rDateTime[i]>"2020-5-26" && bottleDataRaw$rDateTime[i]<"2020-10-15"){bottleDataRaw$precipitation[i] = "wet"}
  else if (bottleDataRaw$rDateTime[i]>"2019-5-26" && bottleDataRaw$rDateTime[i]<"2019-10-15"){bottleDataRaw$precipitation[i] = "wet"}
  else if (bottleDataRaw$rDateTime[i]>"2018-5-26" && bottleDataRaw$rDateTime[i]<"2018-10-15"){bottleDataRaw$precipitation[i] = "wet"}
  else {bottleDataRaw$precipitation[i] = "dry"}
}

############bottle calculations and data grabs from str and seafet##############
bottleIter = nrow(bottleDataRaw);#count number of rows in bottleDataRaw
#look up seafet pH (internal) and str temperature for each bottle, where present
for(i in 1:bottleIter){
    if (bottleDataRaw$site[i] == "em") {#determine site
    phLookupLine = match(bottleDataRaw$rDateTime30Rounded[i], phEm$rDateTime)#identify seafet row number for that site and time
    bottleDataRaw$phInt[i] = phEm$phInt[phLookupLine]#grab the pH from that line and put it into the bottle data table
    bottleDataRaw$seafetDeployment[i] = phEm$collection[phLookupLine]#grab the seafet deployment number from that line and put it into the bottle data table
    tempLookupLine = match(bottleDataRaw$dateTimeEST5Rounded[i], tempEm$dateTimeEST5Rounded)#identify the str row number for that site and time
    bottleDataRaw$tempSTR[i] = tempEm$temperature[tempLookupLine]#grab the str temperature from that line and put it into the bottle data table
  }
   else if (bottleDataRaw$site[i] == "nmac") {#determine site
    phLookupLine = match(bottleDataRaw$rDateTime30Rounded[i], phNmac$rDateTime)#identify seafet row number for that site and time
    bottleDataRaw$phInt[i] = phNmac$phInt[phLookupLine]#grab the pH from that line and put it into the bottle data table
    bottleDataRaw$seafetDeployment[i] = phNmac$collection[phLookupLine]#grab the seafet deployment number from that line and put it into the bottle data table
    tempLookupLine = match(bottleDataRaw$dateTimeEST5Rounded[i], tempNmac$dateTimeEST5Rounded)#identify the str row number for that site and time
    bottleDataRaw$tempSTR[i] = tempNmac$temperature[tempLookupLine]#grab the str temperature from that line and put it into the bottle data table
  }
   else if (bottleDataRaw$site[i] == "smac") {#determine site
    phLookupLine = match(bottleDataRaw$rDateTime30Rounded[i], phSmac$rDateTime)#identify seafet row number for that site and time
    bottleDataRaw$phInt[i] = phSmac$phInt[phLookupLine]#grab the pH from that line and put it into the bottle data table
    bottleDataRaw$seafetDeployment[i] = phSmac$collection[phLookupLine]#grab the seafet deployment number from that line and put it into the bottle data table
    tempLookupLine = match(bottleDataRaw$dateTimeEST5Rounded[i], tempSmac$dateTimeEST5Rounded)#identify the str row number for that site and time
    bottleDataRaw$tempSTR[i] = tempSmac$temperature[tempLookupLine]#grab the str temperature from that line and put it into the bottle data table
  }
   else if (bottleDataRaw$site[i] == "star") {#determine site
    phLookupLine = match(bottleDataRaw$rDateTime30Rounded[i], phStar$rDateTime)#identify seafet row number for that site and time
    bottleDataRaw$phInt[i] = phStar$phInt[phLookupLine]#grab the pH from that line and put it into the bottle data table
    bottleDataRaw$seafetDeployment[i] = phStar$collection[phLookupLine]#grab the seafet deployment number from that line and put it into the bottle data table
    tempLookupLine = match(bottleDataRaw$dateTimeEST5Rounded[i], tempStar$dateTimeEST5Rounded)#identify the str row number for that site and time
    bottleDataRaw$tempSTR[i] = tempStar$temperature[tempLookupLine]#grab the str temperature from that line and put it into the bottle data table
  }
  else{
    bottleDataRaw$phInt[i] = 0
    bottleDataRaw$tempSTR[i] = 0
  }
}
#determine temperature used for the field from STR or SAS if STR doesn't exist
for(i in 1:bottleIter){
  if(is.na(bottleDataRaw$tempSTR[i])&&is.na(bottleDataRaw$tempField[i])){
    bottleDataRaw$tempSource[i] = NA
  }  
  else if(is.na(bottleDataRaw$tempSTR[i])){
    bottleDataRaw$tempUsed[i] = bottleDataRaw$tempField[i]
    bottleDataRaw$tempSource[i] = "SAS"
  } #if tempSTR is absent, use the temperature of the field (usually sas temp)
  else {
    bottleDataRaw$tempUsed[i] = bottleDataRaw$tempSTR[i]
    bottleDataRaw$tempSource[i] = "STR"
  } #otherwise use temperature from the str  
}
#calculate the pH insitu from the spectrophotometric pH, depth, salinity, and the tenmpUsed for each bottle sample. Note alkalinity is arbitrary here
for(i in 1:bottleIter){
  if (is.na(bottleDataRaw$phSpec25[i])||is.na(bottleDataRaw$tempUsed[i])||is.na(bottleDataRaw$depthM[i])||is.na(bottleDataRaw$salinityDensitometer[i])){bottleDataRaw$phSpecField[i]=NA}
  else{bottleDataRaw$phSpecField[i]= pHinsi(bottleDataRaw$phSpec25[i], ALK=2.4e-3, Tinsi=bottleDataRaw$tempUsed[i], Tlab=25, Pinsi=bottleDataRaw$depthM[i], S=bottleDataRaw$salinityDensitometer)
  }
}
#calculate pCO2, omega, and DIC from TA and pH
#The Lueker et al. (2000) constants for K1 and K2, the Perez and Fraga (1987) constant for Kf and the Dickson (1990) constant for Ks are recommended by Dickson et al. (2007). 
for(i in 1:bottleIter){
  if(is.na(bottleDataRaw$tempUsed[i])){output = NA} #if tempSTR is absent, set output to NA
  else { 
    output = carb(flag=15, var1=bottleDataRaw$taCorrected[i]/1000000, var2=bottleDataRaw$dicCorrected[i]/1000000, S=bottleDataRaw$salinityDensitometer[i],T=bottleDataRaw$tempUsed[i], P=bottleDataRaw$depthM[i], k1k2 = "l", kf = "pf", ks = "d", pHscale = "T")
    bottleDataRaw$phCalc[i] = output[6]
  #calculates dic from pH and 
    outputDic = carb(flag=8, var1=bottleDataRaw$phSpecField[i], var2=bottleDataRaw$taCorrected[i]/1000000,  S=bottleDataRaw$salinityDensitometer[i],T=bottleDataRaw$tempUsed[i], P=bottleDataRaw$depthM[i], k1k2 = "l", kf = "pf", ks = "d", pHscale = "T")
    bottleDataRaw$dicCalc[i] = outputDic[16]
    bottleDataRaw$dicCalc[i] = as.numeric(bottleDataRaw$dicCalc[i])
    bottleDataRaw$pco2[i] = outputDic[13]
    bottleDataRaw$omegaArag[i] = outputDic[18]
    }
}
#create dic check
bottleDataRaw$dicCalc<-(as.numeric(bottleDataRaw$dicCalc)*1000000)
bottleDataRaw$dicDiff = (bottleDataRaw$dicCorrected-bottleDataRaw$dicCalc)
#make graphs for salinity normalization of TA and DIC
bottleDataMonitor <- filter(bottleDataRaw,  purpose == "monitor" )
my.formula <- y ~ x
taEndPointGraph <- ggplot(bottleDataMonitor, aes(x = salinityDensitometer, y = taCorrected)) +
  theme_light()+
  geom_point(aes(color = site, shape = precipitation), size = 3)+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method = 'lm', formula = my.formula)+
  stat_cor(label.x = 31, label.y = 2480) +
  stat_regline_equation(label.x = 31, label.y = 2500)
dicEndPointGraph <- ggplot(bottleDataMonitor, aes(x = salinityDensitometer, y = dicCalc)) +
  theme_light()+
  geom_point(aes(color = site, shape = precipitation), size = 3)+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method = 'lm', formula = my.formula)+
  stat_cor(label.x = 31, label.y = 2380) +
  stat_regline_equation(label.x = 31, label.y = 2400)
taDicPointGraph <- ggplot(bottleDataMonitor, aes(x = taCorrected, y = dicCalc)) +
  theme_light()+
  geom_point(aes(color = site, shape = precipitation), size = 3)+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method = 'lm', formula = my.formula)+
  stat_cor(label.x = 2300, label.y = 2380) +
  stat_regline_equation(label.x = 2300, label.y = 2400)
taDicCalcPointGraph <- ggplot(bottleDataMonitor, aes(x = taCorrected, y = dicCalc)) +
  theme_light()+
  geom_point(aes(color = site, shape = precipitation), size = 3)+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method = 'lm', formula = my.formula)+
  stat_cor(label.x = 2300, label.y = 2380) +
  stat_regline_equation(label.x = 2300, label.y = 2400)
dicDiffTimeGraph <- ggplot(bottleDataMonitor, aes(x = rDateTime, y = dicDiff)) +
  theme_light()+
  geom_point(aes(color = site, shape = precipitation), size = 3)+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dic Issues Graph
ggarrange(taDicPointGraph, taDicCalcPointGraph, dicDiffTimeGraph, ncol = 3, nrow = 1)
ggarrange(taEndPointGraph,dicEndPointGraph, ncol = 2, nrow = 1)
#proceed with Friss salinity normalization from end point graph equations
taEndMember = 2800
dicEndMember = 3200
bottleDataMonitor$ta35Friss = 35*((bottleDataMonitor$taCorrected-taEndMember)/bottleDataMonitor$salinityDensitometer)+taEndMember
bottleDataMonitor$dic35Friss = 35*((bottleDataMonitor$dicCalc-dicEndMember)/bottleDataMonitor$salinityDensitometer)+dicEndMember
#graph TA vs DIC
taDicFrissPointGraph <- ggplot(bottleDataMonitor, aes(x = dic35Friss, y = ta35Friss)) +
  theme_light()+
  geom_point(aes(color = site, shape = precipitation), size = 3)+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method = 'lm', formula = my.formula)+
  stat_cor(label.x = 2000, label.y = 2500) +
  stat_regline_equation(label.x = 2000, label.y = 2520)
#create ta dic graph
taDicFrissPointGraph
#create ta and dc vs salinity graphs
ggarrange(taEndPointGraph,dicEndPointGraph, ncol = 2, nrow = 1)



#create bottle tables for all samples
bottleDataMonitor$pco2<-as.numeric(bottleDataMonitor$pco2)
bottleDataMonitor$omegaArag<-as.numeric(bottleDataMonitor$omegaArag)
bottleSummary <- ddply(bottleDataMonitor, c("site","precipitation"), summarise,
                             countSalinityDensitometer = length(dicCalc),
                             meanSalinityDensitometer = round(mean(salinityDensitometer),2),
                             sdSalinityDensitometer   = round(sd(salinityDensitometer),2),
                             meanSpecPhField = round(mean(phSpecField, na.rm=TRUE),3),
                             sdSpecPhField   = round(sd(phSpecField, na.rm=TRUE),3),
                             meanTaCorrected = round(mean(taCorrected, na.rm=TRUE),1),
                             sdTaCorrected   = round(sd(taCorrected, na.rm=TRUE),1),
                             meanDicCalc = round(mean(dicCalc, na.rm=TRUE),1),
                             sdDicCalc   = round(sd(dicCalc, na.rm=TRUE),2),
                             meanTaFriss = round(mean(ta35Friss, na.rm=TRUE),1),
                             sdTaFriss   = round(sd(ta35Friss, na.rm=TRUE),1),
                             meanDicFriss = round(mean(dic35Friss, na.rm=TRUE),1),
                             sdDicFriss   = round(sd(dic35Friss, na.rm=TRUE),1),                            
                             meanPco2 = round(mean(pco2, na.rm=TRUE),0),
                             sdPco2   = round(sd(pco2, na.rm=TRUE),0),
                             meanOmegaArag = round(mean(omegaArag, na.rm=TRUE),2),
                             sdOmegaArag   = round(sd(omegaArag, na.rm=TRUE),2)
)
bottleSummary

#plot TA and DIC as a function of time
dicTimeGraph <- ggplot(bottleDataMonitor, aes(x = time*24, y = dicCalc)) +
  theme_light()+
  geom_point(aes(color = site, shape = precipitation), size = 3)+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+
  scale_x_continuous(name = "time", limits = c(0,24), breaks = c(0,6,12,18,24))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
taTimeGraph <- ggplot(bottleDataMonitor, aes(x = time*24, y = taCorrected)) +
  theme_light()+
  geom_point(aes(color = site, shape = precipitation), size = 3)+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+
  scale_x_continuous(name = "time", limits = c(0,24), breaks = c(0,6,12,18,24))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggarrange(taTimeGraph,dicTimeGraph, ncol = 2, nrow = 1)


#create bottle tables from midday samples (NOTE WE ARE NOT USING THIS BECAUSE IT SEEMS LIKE IT"S NOT STRONGLY RELATED)
bottleDataMidday <- filter(bottleDataMonitor,  midday == "1" )
bottleDataMidday$pco2<-as.numeric(bottleDataMidday$pco2)
bottleDataMidday$omegaArag<-as.numeric(bottleDataMidday$omegaArag)
bottleMiddaySummary <- ddply(bottleDataMidday, c("site","precipitation"), summarise,
                             countSalinityDensitometer = length(dicCalc),
                             meanSalinityDensitometer = round(mean(salinityDensitometer),2),
                             sdSalinityDensitometer   = round(sd(salinityDensitometer),2),
                             meanSpecPhField = round(mean(phSpecField, na.rm=TRUE),3),
                             sdSpecPhField   = round(sd(phSpecField, na.rm=TRUE),3),
                             meanTaCorrected = round(mean(taCorrected, na.rm=TRUE),1),
                             sdTaCorrected   = round(sd(taCorrected, na.rm=TRUE),1),
                             meanDicCalc = round(mean(dicCalc, na.rm=TRUE),1),
                             sdDicCalc   = round(sd(dicCalc, na.rm=TRUE),2),
                             meanTaFriss = round(mean(ta35Friss, na.rm=TRUE),1),
                             sdTaFriss   = round(sd(ta35Friss, na.rm=TRUE),1),
                             meanDicFriss = round(mean(dic35Friss, na.rm=TRUE),1),
                             sdDicFriss   = round(sd(dic35Friss, na.rm=TRUE),1),                            
                             meanPco2 = round(mean(pco2, na.rm=TRUE),0),
                             sdPco2   = round(sd(pco2, na.rm=TRUE),0),
                             meanOmegaArag = round(mean(omegaArag, na.rm=TRUE),2),
                             sdOmegaArag   = round(sd(omegaArag, na.rm=TRUE),2)
)
bottleMiddaySummary














#calculate the particular seafet offset for each bottle sample
for(i in 1:bottleIter){
  if(is.na(bottleDataRaw$phSpecField[i])){#if spec pH field is absent, use the pH calculated from dic and ta
#    bottleDataRaw$offset[i] =1
    bottleDataRaw$offset[i] = (bottleDataRaw$phInt[i] - as.numeric(bottleDataRaw$phCalc[i]))
    bottleDataRaw$phUsedOffset[i] = as.numeric(bottleDataRaw$phCalc[i])
    bottleDataRaw$offsetSource[i] = "calc"
    
  } 
  else {
#    bottleDataRaw$offset[i] =2
    bottleDataRaw$offset[i] = (bottleDataRaw$phInt[i] - bottleDataRaw$phSpecField[i])
    bottleDataRaw$phUsedOffset[i] = bottleDataRaw$phSpecField[i]
    bottleDataRaw$offsetSource[i] = "spec"
  }
}

############graphs for offset analysis####################
phSiteGraph <- ggplot()+
  geom_line(data = phDataRaw, aes(x = rDateTime, y = phInt,  color = site)) +
  geom_point(data = bottleDataRaw, aes(x = rDateTime, y = phUsedOffset, color = site, shape = offsetSource)) +
  scale_shape_manual(values=c(15, 19))+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+  
  theme_classic()+
  scale_y_continuous(name = "pH Int", limits = c(7.15,8.15), labels = scales::number_format(accuracy = 0.1))+
  scale_x_datetime(date_breaks = "1 month", labels = NULL, limits = c(as.POSIXct("2018-12-15"), as.POSIXct("2020-07-15")))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank())+
  ggtitle("ph Internal") + 
  theme(plot.title = element_text(size = 14, face="bold"))+
  geom_line()
phSiteGraph


#Offset table made from offsetDeploymentNotes
#Sites as columns rows as deployment periods

seafetOffset <- cbind(c(-0.0119,-.0151,NA),c(-.0799,NA,-.0118),c(-.0551,-.0407,.0269),c(NA,-.0776,-.0022))

phDataQAQC <- filter(phDataRaw,  good == "y" )#collect all good data 
phIter = nrow(phDataQAQC);#count number of rows in bottleDataRaw
#look up seafet pH (internal) and str temperature for each bottle, where present
for(i in 1:phIter){
  if ((phDataQAQC$site[i] == "em") && (phDataQAQC$collection[i] == "1")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[1,1]}#determine site
  else if ((phDataQAQC$site[i] == "em") && (phDataQAQC$collection[i] == "2")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[2,1]}#determine site
  else if ((phDataQAQC$site[i] == "em") && (phDataQAQC$collection[i] == "3")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[3,1]}#determine site
  else if ((phDataQAQC$site[i] == "nmac") && (phDataQAQC$collection[i] == "1")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[1,2]}#determine site
  else if ((phDataQAQC$site[i] == "nmac") && (phDataQAQC$collection[i] == "2")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[2,2]}#determine site
  else if ((phDataQAQC$site[i] == "nmac") && (phDataQAQC$collection[i] == "3")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[3,2]}#determine site
  else if ((phDataQAQC$site[i] == "smac") && (phDataQAQC$collection[i] == "1")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[1,3]}#determine site
  else if ((phDataQAQC$site[i] == "smac") && (phDataQAQC$collection[i] == "2")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[2,3]}#determine site
  else if ((phDataQAQC$site[i] == "smac") && (phDataQAQC$collection[i] == "3")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[3,3]}#determine site
  else if ((phDataQAQC$site[i] == "star") && (phDataQAQC$collection[i] == "1")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[1,4]}#determine site
  else if ((phDataQAQC$site[i] == "star") && (phDataQAQC$collection[i] == "2")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[2,4]}#determine site
  else if ((phDataQAQC$site[i] == "star") && (phDataQAQC$collection[i] == "3")) {phDataQAQC$phFinal[i]=phDataQAQC$phInt[i]-seafetOffset[3,4]}#determine site
  else{phDataQAQC$phFinal[i]=NA}
}

############graphs for final figure####################
phSiteGraph <- ggplot()+
  geom_line(data = phDataQAQC, aes(x = rDateTime, y = phFinal,  color = site)) +
  geom_point(data = bottleDataRaw, aes(x = rDateTime, y = phUsedOffset, color = site, shape = offsetSource)) +
  scale_shape_manual(values=c(15, 19))+
  scale_color_manual(values=c("#D69EA0","#D2EEEA","#67A9B6","#2A5676"))+  
  theme_classic()+
  scale_y_continuous(name = "pH Int", limits = c(7.5,8.2), breaks = c(7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2), labels = scales::number_format(accuracy = 0.1))+
  scale_x_datetime(date_breaks = "1 month", labels = date_format("%m/%y"), limits = c(as.POSIXct("2018-12-15"), as.POSIXct("2020-07-15")))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank())+
  ggtitle("ph Internal") + 
  theme(plot.title = element_text(size = 14, face="bold"))+
  geom_line()
phSiteGraph


############summary tables####################
seafetSummary <- ddply(phDataQAQC, c("site","precipitation"), summarise,
                       meanPhFinal = round(mean(phFinal),3),
                       sdPhFinal   = round(sd(phFinal),3),
                       minPhFinal = round(min(phFinal),3),
                       maxPhFinal = round(max(phFinal),3))
seafetSummary








###########RUN STATISTICS#######################

