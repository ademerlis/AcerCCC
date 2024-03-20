# playing with openair package
library(RColorBrewer)
library(scales)
library(lubridate)
library(openair)

fileEmerald <- "./2102069_Allyson-2022-08_(0)_Current.csv"
fileMacNorth <- "./2102068_Allyson-2022-08_(0)_Current.csv"

dataEmerald <- read.csv(fileEmerald, header = TRUE)
dataEmerald$date = ymd_hms(paste0(dataEmerald$Date, dataEmerald$Time), tz = "UTC")
dataEmerald <- dataEmerald[, c("date","Speed..cm.s.","Heading..degrees.")]
colnames(dataEmerald) <- c("date","ws","wd")
dataEmerald <- subset(dataEmerald,
                      (date >= ymd_hms("2022/08/11 14:30:00") & date <= ymd_hms("2022/09/30 13:15:00")) |
                        (date >= ymd_hms("2022/10/26 13:30:00") & date <= ymd_hms("2023/11/20 16:45:00"))
)

dataMacNorth <- read.csv(fileMacNorth, header = TRUE)
dataMacNorth$date = ymd_hms(paste0(dataMacNorth$Date, dataMacNorth$Time), tz = "UTC")
dataMacNorth <- dataMacNorth[, c("date","Speed..cm.s.","Heading..degrees.")]
colnames(dataMacNorth) <- c("date","ws","wd")
dataMacNorth <- subset(dataMacNorth,
                       (date >= ymd_hms("2022/08/11 19:00:00") & date <= ymd_hms("2023/06/06 09:45:00")) |
                         (date >= ymd_hms("2023/06/19 08:15:00") & date <= ymd_hms("2023/07/08 12:15:00")) |
                         (date >= ymd_hms("2023/08/01 11:00:00") & date <= ymd_hms("2023/08/12 07:00:00")) |
                         (date >= ymd_hms("2023/09/05 21:45:00") & date <= ymd_hms("2023/11/20 19:15:00"))
)

windRose(dataEmerald, ws.int=3, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="Emerald Reef")
windRose(dataMacNorth, ws.int=6, cols = brewer_pal(palette="Spectral", direction=-1)(6), breaks=6, paddle=FALSE, key.header="(cm/s)", key.footer="MacArthur North")
