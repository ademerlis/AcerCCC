# playing with openair package
library(RColorBrewer)
library(scales)
library(lubridate)
library(openair)

fileRbow <- "~/Documents/aACCRETE/instruments/aRecovery/20201015-tilts-recovery/raw-data/tilt-rainbow-2006057/tilt-rainbow-sn2006057-current.csv"
fileNMac <- "~/Documents/aACCRETE/instruments/aRecovery/20201015-tilts-recovery/raw-data/tilt-northmac-1604207/tilt-northmac-sn1604207-current.csv"
fileCure <- "~/Documents/aACCRETE/instruments/aRecovery/20201015-tilts-recovery/raw-data/tilt-cures-2006058/tilt-cures-sn2006058-current.csv"

dataRbow <- read.csv(fileRbow, header = TRUE)
dataRbow$date = ymd_hms(paste0(dataRbow$Date, dataRbow$Time), tz = "UTC")
dataRbow <- dataRbow[, c("date","Speed..cm.s.","Heading..degrees.")]
colnames(dataRbow) <- c("date","ws","wd")
dataRbow <- subset(dataRbow, date >= ymd_hms("2020/08/05 16:15:00") & date <= ymd_hms("2020/10/13 13:45:00"))

dataNMac <- read.csv(fileNMac, header = TRUE)
dataNMac$date = ymd_hms(paste0(dataNMac$Date, dataNMac$Time), tz = "UTC")
dataNMac <- dataNMac[, c("date","Speed..cm.s.","Heading..degrees.")]
colnames(dataNMac) <- c("date","ws","wd")
dataNMac <- subset(dataNMac, date >= ymd_hms("2020/08/05 13:45:00") & date <= ymd_hms("2020/10/13 14:30:00"))

dataCure <- read.csv(fileCure, header = TRUE)
dataCure$date = ymd_hms(paste0(dataCure$Date, dataCure$Time), tz = "UTC")
dataCure <- dataCure[, c("date","Speed..cm.s.","Heading..degrees.")]
colnames(dataCure) <- c("date","ws","wd")
dataCure <- subset(dataCure, date >= ymd_hms("2020/09/18 15:15:00") & date <= ymd_hms("2020/10/13 15:45:00"))

windRose(dataRbow, ws.int=5, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="Rainbow Reef")
windRose(dataNMac, ws.int=5, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="North MacArthur")
windRose(dataCure, ws.int=5, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="CURES Nursery")
