#load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(zoo)

#snow timing
#precipitation

#load climate data
setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/DataForAnalysis/climate/")
clim= read.csv("AlexanderClimateAll_filled_May2022.csv")

clim= clim[clim$Site=="C1",]
#clim= clim[clim$Site=="NOAA",]

#add precipitation at c1
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nwt.184.5
setwd("/Volumes/GoogleDrive/Shared drives/RoL_FitnessConstraints/projects/Extremes/data/C1_precip/")
pre= read.csv("c1_infilled_precip_daily.tk.data.csv")
date= as.POSIXlt(pre$date, format = "%Y-%m-%d")
pre$Ordinal= date$yday +1
pre$yr_ord= paste(pre$year,pre$Ordinal, sep="_")
#add to clim
clim$yr_ord= paste(clim$Year,clim$Julian, sep="_")
clim$pre= pre$precip[match(clim$yr_ord, pre$yr_ord)]

#load Niwot snowtel
setwd("/Volumes/GoogleDrive/Shared drives/RoL_FitnessConstraints/projects/Extremes/data/NiwotSnowtel/")
snow= read.csv("NiwotSnowtel.csv")
date= as.POSIXlt(snow$Date, format = "%m/%d/%y")
snow$Ordinal= date$yday +1
snow$Year= format(as.Date(snow$Date, format = "%m/%d/%y"),"%Y")
#add to clim
snow$yr_ord= paste(snow$Year,snow$Ordinal, sep="_")
clim$snow= snow$Snow.Water.Equivalent..in..Start.of.Day.Values[match(clim$yr_ord, snow$yr_ord)]

#code season
clim$season<- "summer"
clim$season[clim$Julian <166]<- "spring"
clim$season[clim$Julian >227]<- "fall"

#mean temps by year
clim.ann= aggregate(clim, by=list(clim$Year), FUN="mean", na.rm=TRUE)
ggplot(data=clim.ave, aes(x=Year, y =Max)) + geom_point()
ggplot(data=clim.ave, aes(x=Year, y =pre)) + geom_point()
ggplot(data=clim.ave, aes(x=Year, y =snow)) + geom_point()
#2012 low average snowpack

#seasonal data
clim.ave= aggregate(clim, by=list(clim$Year, clim$season), FUN="mean", na.rm=TRUE)
names(clim.ave)[1:2]<-c("year","season")
ggplot(data=clim.ave[,2:7], aes(x=Year, y =Max,color=season)) + geom_line()
ggplot(data=clim.ave[,c(2:7,11)], aes(x=Year, y =pre,color=season)) + geom_line()
ggplot(data=clim.ave[,c(2:7,11)], aes(x=Year, y =snow,color=season)) + geom_line()
#2012 spring is warm, not other seasons; wet summer
#2012 low spring snow

#add seasonal data 
match1= match(clim.ann$Year, clim.ave[clim.ave$season=="spring","Year"])
clim.ann$SpringPre= clim.ave[clim.ave$season=="spring","pre"][match1]
clim.ann$SpringSnow= clim.ave[clim.ave$season=="spring","snow"][match1]

match1= match(clim.ann$Year, clim.ave[clim.ave$season=="summer","Year"])
clim.ann$SummerPre= clim.ave[clim.ave$season=="summer","pre"][match1]
clim.ann$SummerSnow= clim.ave[clim.ave$season=="summer","snow"][match1]

#-------------------------------
#calculate CLIMDEX indexes
library(climdex.pcic)
library(PCICt)

#estimate climate extremes for each site
years= unique(clim$Year)
sites=c("CHA", "A1", "B1", "C1")

#data storage
clim.ext= matrix(data=NA, nrow=length(years), ncol=length(sites)*2)
clim.ext= as.data.frame(clim.ext)
colnames(clim.ext)=c("CHA_tx90", "CHA_wsdi","A1_tx90", "A1_wsdi","B1_tx90", "B1_wsdi","C1_tx90", "C1_wsdi")

#for sites
for(site in 1:length(sites)){
  clim.sub= clim[which(clim$Site==sites[site] ),]

## Create a climdexInput object from some data already loaded in and
## ready to go.
## Parse the dates into PCICt.
tmax.dates <- as.PCICt(do.call(paste, clim.sub[,c("Year",
                  "Julian")]), format="%Y %j", 
                  cal="gregorian")
tmin.dates <- as.PCICt(do.call(paste, clim.sub[,c("Year",
                                              "Julian")]), format="%Y %j", 
                       cal="gregorian")
prec.dates <- as.PCICt(do.call(paste, clim.sub[,c("Year",
                                              "Julian")]), format="%Y %j", 
                       cal="gregorian")

## Load the data in.
ci <- climdexInput.raw(clim.sub$Max,
                       clim.sub$Min,prec=NULL,
                       tmax.dates, tmin.dates, prec.dates=NULL,
                       base.range=c(1958, 2015),
                       quantiles=NULL, max.missing.days = c(annual = 50,
                                                            monthly = 5) )
#https://search.r-project.org/CRAN/refmans/climdex.pcic/html/climdexInput.raw.html

##Absolute metrics		
#TXx	Annual maximum temperature	°C	Annual maximum of daily TX
#TNn	Annual minimum temperature	°C	Annual minimum of daily TN
#GSL	Growing season length	days	
#DTR	Diurnal temperature range	°C	Mean of daily TX–TN

##Relative metrics
#WSI	Warm spell incidence	%	% Years with ≥6 consecutive days when TX > TX90
#CSI	Cold spell incidence	%	% Years with ≥6 consecutive nights when TN < TN10
#TX90p	Warm days	%	% Days when TX > TX90
#TN10p	Cold nights	%	% Nights when TN < TN10

txx<- climdex.txx(ci, freq = c("monthly"))
tnn<- climdex.tnn(ci, freq = c("monthly"))

tx10p <- climdex.tx10p(ci, freq = c("monthly"))
tx90p <- climdex.tx90p(ci, freq = c("monthly"))

#plot
tx90= as.data.frame(cbind(names(tx90p), tx90p))
tx90$year= format(as.Date(tx90$V1, format="%Y-%M"),"%Y")
tx90$month= as.numeric(substr(tx90$V1, 6, 7))
tx90$tx90p= as.numeric(tx90$tx90p)
tx90$year= as.numeric(tx90$year)

sum.tx90= tx90[tx90$month %in% 6:8,]
sum.tx90= aggregate(sum.tx90, by=list(sum.tx90$year), FUN="mean", na.rm=TRUE)
sum.tx90$year= sum.tx90$Group.1

#ggplot(data=sum.tx90, aes(x=year, y =tx90p )) + geom_point()

#WSDI, Warm spell duration index: Annual count of days with at least 6 consecutive days when TX > 90th percentile

fd <- number.days.op.threshold(ci@data$tmin,
                               ci@date.factors$annual, 0, "<")
fd= as.data.frame(cbind(names(fd), fd))

#store data
match1= match(years, sum.tx90$year)
clim.ext[,(site-1)*2+1]= sum.tx90$tx90p[match1]

match1= match(years, fd$V1)
clim.ext[,(site-1)*2+2]= fd$fd[match1]

} #end loop sites

#write out
#clim.ext$year=years
clim.ann= cbind(clim.ann, clim.ext)

#Save yearly data
setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/DataForAnalysis/climate/")
write.csv(clim.ann, "AlexanderYearlyClimate.csv")

#OTHER METRICS
#ci@quantiles$tmax$inbase$q90
#select.blocks.gt.length(d, n, na.value = FALSE)
#spell.length.max(daily.prec, date.factor, threshold, op, spells.can.span.years)
#threshold.exceedance.duration.index(daily.temp, date.factor, jdays, thresholds,
#                                    op = ">", min.length = 6, spells.can.span.years = # TRUE, max.missing.days)


#-----------------------------


