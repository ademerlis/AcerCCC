# playing with openair package
library(RColorBrewer)
library(scales)
library(lubridate)
library(openair)

fileRbow <- "~/Documents/aACCRETE/instruments/aRecovery/20201015-tilts-recovery/raw-data/tilt-rainbow-2006057/tilt-rainbow-sn2006057-current.csv"
fileNMac <- "~/Documents/aACCRETE/instruments/aRecovery/20201015-tilts-recovery/raw-data/tilt-northmac-1604207/tilt-northmac-sn1604207-current.csv"
fileCure <- "~/Documents/aACCRETE/instruments/aRecovery/20201015-tilts-recovery/raw-data/tilt-cures-2006058/tilt-cures-sn2006058-current.csv"
fileCur2 <- "~/Documents/aACCRETE/instruments/aRecovery/20211209-urban-tilt/1901102_urban-2020-10_(0)_Current.csv"
fileBuoy <- "~/Documents/aACCRETE/instruments/aRecovery/20220127-keys-recv/tiltmeter/2102066_Cheeca-2021-08_(0)_Current.csv"

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

dataCur2 <- read.csv(fileCur2, header = TRUE)
dataCur2$date = ymd_hms(paste0(dataCur2$Date, dataCur2$Time), tz = "UTC")
dataCur2 <- dataCur2[, c("date","Speed..cm.s.","Heading..degrees.")]
colnames(dataCur2) <- c("date","ws","wd")
dataCur2 <- subset(dataCur2, date >= ymd_hms("2020/10/13 16:30:00") & date <= ymd_hms("2021/12/09 18:00:00"))

dataBuoy <- read.csv(fileBuoy, header = TRUE)
dataBuoy$date = ymd_hms(paste0(dataBuoy$Date, dataBuoy$Time), tz = "UTC")
dataBuoy <- dataBuoy[, c("date","Speed..cm.s.","Heading..degrees.")]
colnames(dataBuoy) <- c("date","ws","wd")
dataBuoy <- subset(dataBuoy, date >= ymd_hms("2021/08/31 14:15:00") & date <= ymd_hms("2022/01/27 15:15:00"))

windRose(dataRbow, ws.int=5, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="Rainbow Reef")
windRose(dataNMac, ws.int=5, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="North MacArthur")
windRose(dataCure, ws.int=5, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="CURES Nursery")
windRose(dataCur2, ws.int=5, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="CURES Nursery")
windRose(dataBuoy, ws.int=5, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="Cheeca Buoy")
