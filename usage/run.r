## ======================================================================
## Step 1: Download all log files
## ======================================================================
 
# Here's an easy way to get all the URLs in R
start <- as.Date('2014-04-01')
today <- as.Date('2014-04-28')
 
all_days <- seq(start, today, by = 'day')
 
year <- as.POSIXlt(all_days)$year + 1900
urls <- paste0('http://cran-logs.rstudio.com/', year, '/', all_days, '.csv.gz')
 
# only download the files you don't have:
missing_days <- setdiff(as.character(all_days), tools::file_path_sans_ext(dir("CRANlogs"), TRUE))
 
dir.create("CRANlogs")
for (i in 1:length(missing_days)) {
  print(paste0(i, "/", length(missing_days)))
  download.file(urls[i], paste0('CRANlogs/', missing_days[i], '.csv.gz'))
}

 
 
## ======================================================================
## Step 2: Load single data files into one big data.table
## ======================================================================
 
file_list <- list.files("CRANlogs", full.names=TRUE)
 
logs <- list()
for (file in file_list) {
	print(paste("Reading", file, "..."))
	logs[[file]] <- read.table(file, header = TRUE, sep = ",", quote = "\"",
         dec = ".", fill = TRUE, comment.char = "", as.is=TRUE)
}
 
# rbind together all files
library(data.table)
dat <- rbindlist(logs)
 
# add some keys and define variable types
dat[, date:=as.Date(date)]
dat[, package:=factor(package)]
dat[, country:=factor(country)]
dat[, weekday:=weekdays(date)]
dat[, week:=strftime(as.POSIXlt(date),format="%Y-%W")]
 
setkey(dat, package, date, week, country)
 
save(dat, file="CRANlogs/CRANlogs.RData")
 
# for later analyses: load the saved data.table
# load("CRANlogs/CRANlogs.RData")


## ======================================================================
## Step 3: Analyze it!
## ======================================================================
 
library(ggplot2)
library(plyr)
 
str(dat)
 
# Overall downloads of packages
d1 <- dat[, length(week), by=package]
d1 <- d1[order(V1), ]
d1[package=="PReMiuM", ]
d1[package=="mixcat", ]
 
# plot 1: Compare downloads of selected packages on a weekly basis
agg1 <- dat[J(c("PReMiuM", "mixcat")), length(unique(ip_id)), by=c("week", "package")]
 
ggplot(agg1, aes(x=week, y=V1, color=package, group=package)) + geom_line() + ylab("Downloads") + theme_bw() + theme(axis.text.x  = element_text(angle=90, size=8, vjust=0.5))
 
 
agg1 <- dat[J(c("PReMiuM", "mixcat", "mixcat")), length(unique(ip_id)), by=c("week", "package")]
 
ggplot(agg1, aes(x=week, y=V1, color=package, group=package)) + geom_line() + ylab("Downloads") + theme_bw() + theme(axis.text.x  = element_text(angle=90, size=8, vjust=0.5))
