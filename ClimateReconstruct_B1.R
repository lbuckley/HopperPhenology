#load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(zoo)
library(reshape)

#Reconstruct A1 data for 2009, 2010

#--------------------------------------
fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled.csv")

#Use NOAA, C1, Gross reservoir
#Don't use reconstructed B1 data

#load Gross Reservoir data
setwd( paste(fdir, "climate\\Recent\\", sep="") )   
clim.gr= read.csv("COOP53629_2005_2011.csv")

#Convert F to C
clim.gr$TmaxC= (clim.gr$TMAX-32)*5/9
clim.gr$TminC= (clim.gr$TMIN-32)*5/9

#calculate Julian
tmp <- as.POSIXlt(as.character(clim.gr$DATE), format = "%Y%m%d")
clim.gr$Julian= tmp$yday+1
clim.gr$Year= as.numeric(substr(as.character(clim.gr$DATE), 1,4))
clim.gr$Site="GR"
clim.gr$sjy= paste(clim.gr$Site, clim.gr$Julian, clim.gr$Year, sep="_")

#set up data frame with all combinations
clim1.gr = expand.grid(Site="GR", Julian = 60:243, Year = 2005:2011)
#set up columns
clim1.gr$Max=NA; clim1.gr$Min=NA; clim1.gr$Mean=NA
#columns for matching
clim1.gr$sjy= paste(clim1.gr$Site, clim1.gr$Julian, clim1.gr$Year, sep="_")

#match
match1= match(clim1.gr$sjy, clim.gr$sjy)
matched= which(!is.na(match1))

clim1.gr[matched,c("Max","Min")]= clim.gr[match1[matched], c("TmaxC","TminC")]

#add
clim= rbind(clim[,2:8],clim1.gr)

#==========================================
#PLOT RELATIONSHIPS

clim1= clim[clim$Year %in%2005:2008 & clim$Site %in% c("A1","C1","NOAA","GR"),]

#reformat to wide format for Julian, Year by Site 
clim.min= cast(clim1, Julian +Year ~ Site, value="Min")
clim.max= cast(clim1, Julian +Year ~ Site, value="Max")

#MIN
#plots
clim2= clim.min
clim2= clim.max
ggplot(data=clim2, aes(x=A1, y = NOAA, color=Year))+geom_point() +theme_bw()
ggplot(data=clim2, aes(x=A1, y = C1, color=Year))+geom_point() +theme_bw()
ggplot(data=clim2, aes(x=A1, y = GR, color=Year))+geom_point() +theme_bw()

#models
clim2= clim.min
min.noaa= summary(lm(A1~NOAA, data=clim2))
min.c1= summary(lm(A1~C1, data=clim2))
min.gr= summary(lm(A1~GR, data=clim2)) #SKIP WORST FIT
clim2= clim.max
max.noaa= summary(lm(A1~NOAA, data=clim2))
max.c1= summary(lm(A1~C1, data=clim2))
max.gr= summary(lm(A1~GR, data=clim2)) #SKIP WORST FIT

#------------------
#Predict 2009 and 2010
clim1= clim[clim$Year %in%2009:2010 & clim$Site %in% c("A1","C1","NOAA","GR"),]

#reformat to wide format for Julian, Year by Site 
clim.min= cast(clim1, Julian +Year ~ Site, value="Min")
clim.max= cast(clim1, Julian +Year ~ Site, value="Max")

#predict min
clim.min$A1.pnoaa= min.noaa$coefficients[1,1] +min.noaa$coefficients[2,1]*clim.min$NOAA
clim.min$A1.pc1= min.c1$coefficients[1,1] +min.c1$coefficients[2,1]*clim.min$C1
#weight by r2
clim.min$A1=(min.noaa$adj.r.squared*clim.min$A1.pnoaa + min.c1$adj.r.squared*clim.min$A1.pc1)/(min.noaa$adj.r.squared+min.c1$adj.r.squared)

#predict max
clim.max$A1.pnoaa= max.noaa$coefficients[1,1] +max.noaa$coefficients[2,1]*clim.max$NOAA
clim.max$A1.pc1= max.c1$coefficients[1,1] +max.c1$coefficients[2,1]*clim.max$C1
#weight by r2
clim.max$A1=(max.noaa$adj.r.squared*clim.max$A1.pnoaa + max.c1$adj.r.squared*clim.max$A1.pc1)/(max.noaa$adj.r.squared+max.c1$adj.r.squared)

#------------------
#write out reconstructed data

#format
clim.min$jy= paste(clim.min$Julian, clim.min$Year, sep="_")
clim.max$jy= paste(clim.max$Julian, clim.max$Year, sep="_")

#match
match1= match(clim.min$jy, clim.max$jy)
matched= which(!is.na(match1))

clim.min$Max=NA
clim.min[matched,"A1.max"]= clim.max[match1[matched], "A1"]
clim2= clim.min[,c("Julian","Year","A1","A1.max")]
names(clim2)[3:4]=c("Min","Max")

setwd( paste(fdir, "climate\\Recent\\", sep="") )   
write.csv(clim2, "A1_2009_2010_reconstruct.csv")
