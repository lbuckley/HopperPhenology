library(ggplot2)
library(plyr)
library(dplyr)

#READ DATA
setwd("C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\TPCfield\\data")
#setwd("C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\TPCfield\\data")

hop= read.csv("HoppingData.csv")
feed=  read.csv("FeedingData.csv")
pbt= read.csv("datall_PBT.csv")

sites= c("Redfox", "A1", "B1", "C1", "D1")  
elevs= c(1574, 2195, 2591, 3048, 3739)
specs=c("M. dodgei", "C. pellucida", "M. sanguinipes", "A. clavatus")

count= function(x) length(na.omit(x))
se= function(x) sd(x)#/sqrt(length(x))

#------------
#HOPPING
hop= hop[hop$Species=="dodgei",]

#CONVERT FROM FT TO M
hop$dist =hop$dist*0.3048
hop$elev= as.factor(hop$elev)

#group
hop2= hop %>% group_by(temp,elev,Sex) %>% summarise(dist=mean(dist), Mass=mean(Mass), FemurLength=mean(FemurLength), dist.se=se(dist) )

#PLOT
pd <- position_dodge(0.1)
ggplot(data=hop2, aes(x=temp, y = dist,color=elev, shape=Sex))+geom_point()+geom_line()+ 
  geom_errorbar(aes(ymin=dist-dist.se, ymax=dist+dist.se), width=.1, position=pd)

#FOR FITTING PERFORMANCE CURVES
TPC.gaus= function(T, b, c){1/c*exp(-0.5*((T-b)/c)^2)}
TPC.gausgomp= function(T, To, rho=0.9, sigma, Fmax) Fmax*exp(-exp(rho*(T-To)-6)-sigma*(T-To)^2) # rho and sigma determine, the thermal sensitivity of feeding at temperatures above and below Topt, respectively

#---------------------
#FEEDING
feed= feed[feed$Species=="dodgei",]

#Add elevation
feed$Elev=elevs[match(feed$Site,sites)]

#Adjust for mass and normalize pixel #
feed$Area_8hr_nomass= feed$Area_8hr
feed$Area_2hr_nomass= feed$Area_2hr
feed$Area_8hr= feed$Area_8hr/feed$Mass /10^6
feed$Area_2hr= feed$Area_2hr/feed$Mass /10^6

#group
feed2= feed %>% group_by(Temp,Elev,Sex) %>% summarise(Area_2hr=mean(Area_2hr),Area_8hr=mean(Area_8hr),Area_2hr_nomass=mean(Area_2hr_nomass),Area_8hr_nomass=mean(Area_8hr_nomass), Mass=mean(Mass), FemurLength=mean(FemurLength) )
feed2$Elev= as.factor(feed2$Elev)

#PLOT
ggplot(data=feed2, aes(x=Temp, y = Area_8hr_nomass,color=Elev, shape=Sex))+geom_point()+geom_line()

#---------------------
#PBT
pbt= pbt[pbt$Spec=="M. dodgei",]

#group
pbt2= pbt %>% group_by(elev,sex) %>% summarise(PBT=mean(PBT, na.rm=TRUE))
pbt$elev= as.factor(pbt$elev)

#PLOT
ggplot(data=pbt2, aes(x=elev, y = PBT,color=elev, shape=sex))+geom_point()+geom_line()



