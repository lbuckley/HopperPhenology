library(ggplot2)
library(plyr)
library(dplyr)
library(grid)

fdir= "C:\\Users\\Buckley\\Google Drive\\Work\\GrasshopperGenetics\\"

#READ DATA
#setwd("C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\TPCfield\\data")
setwd("C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\TPCfield\\data")

hop= read.csv("HoppingData.csv")
feed=  read.csv("FeedingData.csv")
pbt= read.csv("datall_PBT.csv")
egg= read.csv("LevyNufioData.csv")

#change dodgei to boulderensis
hop$Species= gsub("dodgei", "boulderensis", hop$Species)
feed$Species= gsub("dodgei", "boulderensis", feed$Species)
pbt$Spec= gsub("dodgei", "boulderensis", pbt$Spec)

sites= c("Redfox", "A1", "B1", "C1", "D1")  
elevs= c(1574, 2195, 2591, 3048, 3739)
specs=c("M. boulderensis", "C. pellucida", "M. sanguinipes", "A. clavatus")

count= function(x) length(na.omit(x))
se= function(x) sd(x)/sqrt(length(x))

#------------
#HOPPING
#hop= hop[hop$Species=="dodgei",]

#CONVERT FROM FT TO M
hop$dist =hop$dist*0.3048
hop$elev= as.factor(hop$elev)

#group
hop2= hop %>% group_by(temp,elev,Species) %>% summarise(N= length(dist), dist=mean(dist), Mass=mean(Mass), FemurLength=mean(FemurLength), dist.sd=sd(dist), dist.se = dist.sd / sqrt(N))
#average over sex for now

#PLOT
pd <- position_dodge(0.1)
hop.plot=ggplot(data=hop2, aes(x=temp, y = dist,color=elev))+geom_point()+geom_line()+ 
  geom_errorbar(aes(ymin=dist-dist.se, ymax=dist+dist.se), width=.1, position=pd)+facet_wrap(~Species)+theme_bw()
#, shape=Sex

#FOR FITTING PERFORMANCE CURVES
TPC.gaus= function(T, b, c){1/c*exp(-0.5*((T-b)/c)^2)}
TPC.gausgomp= function(T, To, rho=0.9, sigma, Fmax) Fmax*exp(-exp(rho*(T-To)-6)-sigma*(T-To)^2) # rho and sigma determine, the thermal sensitivity of feeding at temperatures above and below Topt, respectively

#PLOT
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperGenetics\\figures\\")

#min ordinal date by dd, #plot regression
pdf("hop_plot.pdf", height = 5, width = 10)
hop.plot
dev.off()

#---------------------
#FEEDING

#Add elevation
feed$Elev=elevs[match(feed$Site,sites)]

#Adjust for mass and normalize pixel #
feed$Area_8hr_nomass= feed$Area_8hr
feed$Area_2hr_nomass= feed$Area_2hr
feed$Area_8hr= feed$Area_8hr/feed$Mass /10^6
feed$Area_2hr= feed$Area_2hr/feed$Mass /10^6

#group
feed2= feed %>% group_by(Temp,Elev, Species) %>% summarise(Area_2hr=mean(Area_2hr),Area_8hr=mean(Area_8hr),Area_2hr_nomass=mean(Area_2hr_nomass),Area_8hr_nomass=mean(Area_8hr_nomass), Area_8hr_nomass.sd=sd(Area_8hr_nomass), Mass=mean(Mass), FemurLength=mean(FemurLength) )
feed2$Elev= as.factor(feed2$Elev)
#group by Sex for now

#PLOT
ggplot(data=feed2, aes(x=Temp, y = Area_8hr_nomass,color=Elev))+geom_point()+geom_line()+facet_wrap(~Species)+theme_bw()

#---------------------
#PBT

#group
pbt2= pbt %>% group_by(elev,Spec) %>% summarise(PBT=mean(PBT, na.rm=TRUE), PBT.N= length(na.omit(PBT)),PBT.sd=sd(PBT, na.rm=TRUE), PBT.se=PBT.sd/sqrt(PBT.N), Ctmin=mean(Ctmin, na.rm=TRUE), Ctmin.N= length(na.omit(Ctmin)),Ctmin.sd=sd(Ctmin, na.rm=TRUE),Ctmin.se=Ctmin.sd/sqrt(Ctmin.N), Ctmax=mean(Ctmax, na.rm=TRUE), Ctmax.N= length(na.omit(Ctmax)),Ctmax.sd=sd(Ctmax, na.rm=TRUE), Ctmax.se=Ctmax.sd/sqrt(Ctmax.N), mass=mean(mass_g, na.rm=TRUE), mass.N= length(na.omit(mass_g)),mass.sd=sd(mass_g, na.rm=TRUE), mass.se=mass.sd/sqrt(mass.N))

#PLOT
ctmax.plot= ggplot(data=pbt2, aes(x=elev, y = Ctmax, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

pbt.plot= ggplot(data=pbt2, aes(x=elev, y = PBT, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

ctmin.plot= ggplot(data=pbt2, aes(x=elev, y = Ctmin, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="bottom")

mass.plot= ggplot(data=pbt2, aes(x=elev, y = mass, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="bottom")+ 
  geom_errorbar(aes(ymin=mass-mass.se, ymax=mass+mass.se), width=.1, position=pd)

#PLOT
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperGenetics\\figures\\")

pdf("thermtol_plot.pdf", height = 6, width = 8)
grid::grid.draw(rbind(ggplotGrob(ctmax.plot), ggplotGrob(pbt.plot), ggplotGrob(ctmin.plot), size = "last"))
dev.off()

#STATS
specs=c("M. boulderensis", "C. pellucida", "M. sanguinipes", "A. clavatus")
pbt3= pbt[which(pbt$Spec %in% specs[4]),]

mod.mass= lm(pbt3$mass ~ pbt3$elev) 
summary(mod.mass)
mod.ctmax= lm(pbt3$Ctmax ~ pbt3$elev) 
summary(mod.ctmax)
mod.PBT= lm(pbt3$PBT ~ pbt3$elev) 
summary(mod.PBT)
mod.ctmin= lm(pbt3$Ctmin ~ pbt3$elev) 
summary(mod.ctmin)

#CLINES
#boulderensis: mass decrease
#pellucida: CTmin decrease and tendency for PBT to decrease
#sanguinipes: tendency for PBT to decrease
#clavatus: mass decrease, tendency for Ctmax to decrease, PBT increase

#-----------------------------
#Egg plot

#Add elevation
sites1= c("RF", "A1", "B1", "C1", "D1")  
elevs= c(1574, 2195, 2591, 3048, 3739)
egg$elev=elevs[match(egg$Site,sites1)]

ov.plot=ggplot(data=egg, aes(x=elev, y = Novarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=Novarioles-Novarioles.se, ymax=Novarioles+Novarioles.se), width=.1)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
                
clutch.plot=ggplot(data=egg, aes(x=elev, y = clutch,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=clutch-clutch.se, ymax=clutch+clutch.se), width=.1)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

position_dodge(0.1) # move them .05 to the left and right
egg.plot= ggplot(data=egg, aes(x=elev, y = egg.mass,color=Species))+geom_point(position=pd)+geom_line(position=pd)+theme_bw()+ 
  geom_errorbar(aes(ymin=egg.mass-egg.mass.se, ymax=egg.mass+egg.mass.se), width=.1, position=pd)+ theme(legend.position="bottom")

#PLOT
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperGenetics\\figures\\")

pdf("egg_plot.pdf", height = 6, width = 8)
grid::grid.draw(rbind(ggplotGrob(ov.plot), ggplotGrob(clutch.plot), ggplotGrob(egg.plot), size = "last"))
dev.off()

