library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)

fdir= "C:\\Users\\Buckley\\Google Drive\\Work\\GrasshopperGenetics\\"

#READ DATA
#setwd("C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\TPCfield\\data")
setwd("C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\TPCfield\\data")

hop= read.csv("HoppingData.csv")
feed=  read.csv("FeedingData.csv")
pbt= read.csv("datall_PBT.csv")
egg= read.csv("LevyNufioData.csv")
size= read.csv("Levy_Male&FemaleBodySizesGradient.csv")
repro= read.csv("Levy_FemaleGradientDataGrasshopper.csv")

#change dodgei to boulderensis
hop$Species= gsub("dodgei", "boulderensis", hop$Species)
feed$Species= gsub("dodgei", "boulderensis", feed$Species)
pbt$Spec= gsub("dodgei", "boulderensis", pbt$Spec)
repro$Species= gsub("dodgei", "boulderensis", repro$Species)

#change to abbrev species names
size$Species= gsub("boulderensis", "M. boulderensis", size$Species)
size$Species= gsub("clavatus", "A. clavatus", size$Species)
size$Species= gsub("pellucida", "C. pellucida", size$Species)
size$Species= gsub("sanguinipes", "M. sanguinipes", size$Species)

repro$Species= gsub("boulderensis", "M. boulderensis", repro$Species)
repro$Species= gsub("clavatus", "A. clavatus", repro$Species)
repro$Species= gsub("pellucida", "C. pellucida", repro$Species)
repro$Species= gsub("sanguinipes", "M. sanguinipes", repro$Species)

hop$Species= gsub("boulderensis", "M. boulderensis", hop$Species)
hop$Species= gsub("clavatus", "A. clavatus", hop$Species)
hop$Species= gsub("pellucida", "C. pellucida", hop$Species)
hop$Species= gsub("sanguinipes", "M. sanguinipes", hop$Species)

#return to factors
hop$Species= as.factor(hop$Species)
feed$Species= as.factor(feed$Species)
pbt$Spec= as.factor(pbt$Spec) 
repro$Species= as.factor(repro$Species) 
size$Species= as.factor(size$Species) 

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
  geom_errorbar(aes(ymin=dist-dist.se, ymax=dist+dist.se), width=.1, position=pd)+facet_wrap(~Species)+theme_bw()+xlab("Temperature (°C)")+ylab("Hop distance (m)")+ labs(color="Elevation (m)") 
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
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+ylab("CTmax (°C)")

pbt.plot= ggplot(data=pbt2, aes(x=elev, y = PBT, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

ctmin.plot= ggplot(data=pbt2, aes(x=elev, y = Ctmin, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="bottom")+xlab("Elevation (m)")+ylab("CTmin (°C)")+ labs(color="Species") 

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
#Size plot

#group
size2= size %>% group_by(Site,Elevation, Species) %>% summarise(length=mean(Mean.Femur.Length), length.N= length(na.omit(Mean.Femur.Length)),length.sd=sd(Mean.Femur.Length, na.rm=TRUE), length.se=length.sd/sqrt(length.N))
#group by Sex for now

size.plot= ggplot(data=size2, aes(x=Elevation, y = length, color=Species))+geom_point()+geom_line()+theme_bw()+ 
  geom_errorbar(aes(ymin=length-length.se, ymax=length+length.se), width=.1)+ theme(legend.position="none")+ylab("Femur length (mm)")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

#STATS
specs=c("M. boulderensis", "C. pellucida", "M. sanguinipes", "A. clavatus")
size3= size[which(size$Species %in% specs[4]),]

mod.length= lm(size3$Mean.Femur.Length ~ size3$Elevation) 
summary(mod.length)

#-----------------------------
#Repro plot

#group
repro2= repro %>% group_by(Site,Elevation, Species) %>% summarise(clutch=mean(Clutch.Size, na.rm=TRUE), clutch.N= length(na.omit(Clutch.Size)),clutch.sd=sd(Clutch.Size, na.rm=TRUE), clutch.se=clutch.sd/sqrt(clutch.N),
                                                                  egg.mass=mean(Mean.Egg.Mass, na.rm=TRUE), egg.mass.N= length(na.omit(Mean.Egg.Mass)),egg.mass.sd=sd(Mean.Egg.Mass, na.rm=TRUE), egg.mass.se=egg.mass.sd/sqrt(egg.mass.N),
                                                                  Novarioles=mean(Number.Ovarioles, na.rm=TRUE), Novarioles.N= length(na.omit(Number.Ovarioles)),Novarioles.sd=sd(Number.Ovarioles, na.rm=TRUE), Novarioles.se=Novarioles.sd/sqrt(Novarioles.N),
                                                                  NFovarioles=mean(Number.Functioing.Ovarioles, na.rm=TRUE), NFovarioles.N= length(na.omit(Number.Functioing.Ovarioles)),NFovarioles.sd=sd(Number.Functioing.Ovarioles, na.rm=TRUE), NFovarioles.se=NFovarioles.sd/sqrt(NFovarioles.N),
                                                                  PFovarioles=mean(Proportion.Functioning.Ovarioles, na.rm=TRUE), PFovarioles.N= length(na.omit(Proportion.Functioning.Ovarioles)),PFovarioles.sd=sd(Proportion.Functioning.Ovarioles, na.rm=TRUE), PFovarioles.se=PFovarioles.sd/sqrt(PFovarioles.N), 
                                                                  clutch.mass=mean(clutch.weight.g, na.rm=TRUE), clutch.mass.N= length(na.omit(clutch.weight.g)),clutch.mass.sd=sd(clutch.weight.g, na.rm=TRUE), clutch.mass.se=clutch.mass.sd/sqrt(clutch.mass.N) )

#--------------
ov.plot=ggplot(data=repro2, aes(x=Elevation, y = Novarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=Novarioles-Novarioles.se, ymax=Novarioles+Novarioles.se), width=.1)

clutch.plot=ggplot(data=repro2, aes(x=Elevation, y = clutch,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=clutch-clutch.se, ymax=clutch+clutch.se), width=.1)

egg.mass.plot=ggplot(data=repro2, aes(x=Elevation, y = egg.mass*1000,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=egg.mass-egg.mass.se, ymax=egg.mass+egg.mass.se), width=.1)+ylab("Egg mass (mg)")+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

clutch.mass.plot=ggplot(data=repro2, aes(x=Elevation, y = clutch.mass,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=clutch.mass-clutch.mass.se, ymax=clutch.mass+clutch.mass.se), width=.1)+ylab("Clutch mass (g)")+xlab("Elevation (m)")

Fov.plot=ggplot(data=repro2, aes(x=Elevation, y = NFovarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=NFovarioles-NFovarioles.se, ymax=NFovarioles+NFovarioles.se), width=.1)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+ylab("N functional ovarioles")

PFov.plot=ggplot(data=repro2, aes(x=Elevation, y = PFovarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=PFovarioles-PFovarioles.se, ymax=PFovarioles+PFovarioles.se), width=.1)

#------------------
#PLOT
grid::grid.draw(rbind(ggplotGrob(ov.plot), ggplotGrob(Fov.plot),ggplotGrob(PFov.plot),ggplotGrob(clutch.plot), ggplotGrob(clutch.mass.plot),ggplotGrob(egg.mass.plot), size = "last"))

#------------------
#STATS
specs=c("M. boulderensis", "C. pellucida", "M. sanguinipes", "A. clavatus")
repro3= repro[which(repro$Species %in% specs[4]),]

mod.ov= lm(repro3$Number.Ovarioles ~ repro3$Elevation) 
summary(mod.ov)

mod.fov= lm(repro3$Number.Functioing.Ovarioles ~ repro3$Elevation) 
summary(mod.fov)

mod.pov= lm(repro3$Proportion.Functioning.Ovarioles ~ repro3$Elevation) 
summary(mod.pov)

mod.clutch= lm(repro3$Clutch.Size ~ repro3$Elevation) 
summary(mod.clutch)

mod.egg.mass= lm(repro3$Mean.Egg.Mass ~ repro3$Elevation) 
summary(mod.egg.mass)

mod.clutch.mass= lm(repro3$clutch.weight.g ~ repro3$Elevation) 
summary(mod.clutch.mass)

#=====================================
#PHENOTYPE FIG

#Length
#CTmax
#CTmin

#FOV
#Egg.mass
#Clutch.mass

#TPC


lay <- rbind(c(1,4),
             c(2,5),
             c(3,6),
             c(7,7),
             c(8,8))

#get legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(ctmin.plot)

setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperGenetics\\figures\\")
pdf("phenotype_plot.pdf", height = 10, width = 8)
grid.arrange(size.plot, ctmax.plot, ctmin.plot+ theme(legend.position = 'none'), Fov.plot, egg.mass.plot, clutch.mass.plot, legend, hop.plot, layout_matrix = lay, widths = c(1,1), 
             heights=c(0.5,0.5,0.5,0.1,1))
dev.off()
