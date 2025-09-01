# IMPORT LIBRARIES
library(dlnm)
library(epiDisplay); library(lubridate) ; library(dplyr); library(tidyr); library(plyr); library(scales)
library(splines); library(gnm); library(Epi); library(tsModel);  library(sf) 
library(data.table); library(viridis); library(ggpubr)
library(mgcv) ; library(corrplot); library(mixmeta) ; library(metafor)
library(ggplot2) ; library(patchwork); library(gridExtra)

# Load file with health and exposure time series data 
data <- read.csv('dist_level_data.csv') 
data <- as.data.table(data)

# Define the public holidays in Milan # holidays on weekends removed https://www.timeanddate.com/holidays/italy/2019?hol=1 
milan_holidays <- c(
  "2016-01-01", "2016-01-06", "2016-03-28", "2016-04-25", "2016-06-02", "2016-08-15", "2016-11-01", "2016-12-08", "2016-12-26",
  "2017-01-06", "2017-04-17", "2017-04-25", "2017-05-01", "2017-06-02", "2017-08-15", "2017-11-01", "2017-12-08", "2017-12-25", "2017-12-26",
  "2018-01-01", "2018-04-02", "2018-04-25", "2018-05-01", "2018-08-15", "2018-11-01", "2018-12-25", "2018-12-26",
  "2019-01-01", "2019-04-22", "2019-04-25", "2019-05-01", "2019-08-15", "2019-11-01", "2019-12-25", "2019-12-26"
)

# chose data only till 2019 
data$Date <- as.Date(data$Date)
data <- data[data$Date <= as.Date("2019-12-31"), ]
# Edit date format and add additional time related parameters
data$Date <- as.Date(data$Date)
data$year <- as.factor(substr(data$Date,1,4))
data$month <- as.factor(months(data$Date,abbr=TRUE))
data$day <- format(data$Date, format = "%d")
data$dow <- wday(data$Date) 
data$doy <- yday(data$Date)
# Update the "adult" column
data[, adult := adult + young]
# Generate holiday dates for Italy from 2016 to 2019
data$holidays <- as.integer(data$Date %in% as.Date(milan_holidays))
# define stratum
data[, stratum:=factor(paste(DISTRICT, year, month, dow, sep=":"))] 
data[, RH_03 := frollmean(RH, n = 4, align = "right"), by = DISTRICT] 

# CORRELATION ANALYSIS
setnames(data, "all_CA", "OHCA")
correlation_table <- cor(data[, .(PM25, PM10, NO2, O3, SO2, CO, Temp, RH, OHCA)], method = "spearman")

# Print the correlation table
corrplot(correlation_table, type = 'lower', method = 'square', tl.col = "black", tl.srt = 90, diag = FALSE, 
         cl.pos = 'n', tl.cex=0.8, addCoef.col = "black",  number.digits = 2,
         number.cex = 0.8, col = rev(corrplot::COL2('RdBu')))
