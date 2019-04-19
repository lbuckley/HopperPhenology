library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)

#fdir= "C:\\Users\\Buckley\\Google Drive\\Work\\GrasshopperGenetics\\"
fdir= "/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperGenetics/"

library(RColorBrewer)
bg<-brewer.pal(5, "YlGnBu")
sp.col<-c("yellow",  "#fecc5c", "#fd8d3c","#e31a1c")  

#READ DATA
#setwd("C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\TPCfield\\data")
#setwd("C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\TPCfield\\data")
setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/TPCfield/data")

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
hop$Species= factor(hop$Species, levels=c("A. clavatus","M. boulderensis","C. pellucida", "M. sanguinipes") )
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
hop2= hop %>% group_by(temp,elev,Species) %>% summarise(N= length(dist), dist.m=mean(dist), Mass.m=mean(Mass), FemurLength=mean(FemurLength), dist.sd=sd(dist), dist.se = dist.sd / sqrt(N))
#average over sex for now

#PLOT
pd <- position_dodge(0.1)
hop.plot=ggplot(data=hop2, aes(x=temp, y = dist.m,color=elev, fill=elev))+geom_point(colour="black",pch=21, size=2,show.legend = FALSE)+geom_line(show.legend = FALSE)+ 
  geom_errorbar(aes(ymin=dist.m-dist.se, ymax=dist.m+dist.se), width=.1, position=pd)+facet_wrap(~Species)+theme_bw()+xlab("Temperature (°C)")+ylab("Hop distance (m)")+ labs(color="Elevation (m)")+
  scale_colour_manual(values = bg) +
  scale_fill_manual(values = bg) 
#, shape=Sex

#FOR FITTING PERFORMANCE CURVES
#TPC.gaus= function(T, b, c){1/c*exp(-0.5*((T-b)/c)^2)}
#TPC.gausgomp= function(T, To, rho=0.9, sigma, Fmax) Fmax*exp(-exp(rho*(T-To)-6)-sigma*(T-To)^2) # rho and sigma determine, the thermal sensitivity of feeding at temperatures above and below Topt, respectively

#Just melanoplus and transplant elevations
#hop.plot=ggplot(data=hop2[which(hop2$Species %in% c("M. boulderensis","M. sanguinipes") &  hop2$elev %in% c(2195,2591,3048)),], aes(x=temp, y = dist.m,color=elev))+geom_point()+geom_line()+ 
#  geom_errorbar(aes(ymin=dist.m-dist.se, ymax=dist.m+dist.se), width=.1, position=pd)+facet_wrap(~Species)+theme_bw()+xlab("Temperature (°C)")+ylab("Hop distance (m)")+ labs(color="Elevation (m)") 

#write out data
#write.csv(hop2, "RoLHopdata.csv")

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
pbt2= pbt %>% group_by(elev,Spec) %>% summarise(PBT.m=mean(PBT, na.rm=TRUE), PBT.N= length(na.omit(PBT)),PBT.sd=sd(PBT, na.rm=TRUE), PBT.se=PBT.sd/sqrt(PBT.N), Ctmin.m=mean(Ctmin, na.rm=TRUE), Ctmin.N= length(na.omit(Ctmin)),Ctmin.sd=sd(Ctmin, na.rm=TRUE),Ctmin.se=Ctmin.sd/sqrt(Ctmin.N), Ctmax.m=mean(Ctmax, na.rm=TRUE), Ctmax.N= length(na.omit(Ctmax)),Ctmax.sd=sd(Ctmax, na.rm=TRUE), Ctmax.se=Ctmax.sd/sqrt(Ctmax.N), mass.m=mean(mass_g, na.rm=TRUE), mass.N= length(na.omit(mass_g)),mass.sd=sd(mass_g, na.rm=TRUE), mass.se=mass.sd/sqrt(mass.N))

sp.col<-c("A. clavatus"="yellow",  "M. boulderensis"="#fecc5c", "C. pellucida"="#fd8d3c","M. sanguinipes"="#e31a1c")  
pbt2$Spec= factor(pbt2$Spec, levels=c("A. clavatus","M. boulderensis","C. pellucida","M. sanguinipes"), ordered=TRUE)

#PLOT
ctmax.plot= ggplot(data=pbt2, aes(x=elev, y = Ctmax.m, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+ylab("CTmax (°C)")+ 
  geom_errorbar(aes(ymin=Ctmax.m-Ctmax.se, ymax=Ctmax.m+Ctmax.se), width=.1, position=pd) +xlim(1550,3530)+
  scale_colour_manual(values = sp.col)

pbt.plot= ggplot(data=pbt2, aes(x=elev, y = PBT.m, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+ 
  geom_errorbar(aes(ymin=PBT.m-PBT.se, ymax=PBT.m+PBT.se), width=.1, position=pd)

ctmin.plot= ggplot(data=pbt2, aes(x=elev, y = Ctmin.m, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="bottom")+xlab("Elevation (m)")+ylab("CTmin (°C)")+ labs(color="Species") + 
  geom_errorbar(aes(ymin=Ctmin.m-Ctmin.se, ymax=Ctmin.m+Ctmin.se), width=.1, position=pd) +xlim(1550,3530)+
  scale_colour_manual(values = sp.col)

mass.plot= ggplot(data=pbt2, aes(x=elev, y = mass.m, color=Spec))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="bottom")+ 
  geom_errorbar(aes(ymin=mass.m-mass.se, ymax=mass.m+mass.se), width=.1, position=pd)

#PLOT
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperGenetics/figures/")
#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperGenetics\\figures\\")
grid::grid.draw(rbind(ggplotGrob(ctmax.plot), ggplotGrob(pbt.plot), ggplotGrob(ctmin.plot), size = "last"))

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
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +xlim(1550,3530)+
  scale_colour_manual(values = sp.col)

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
sp.col<-c("A. clavatus"="yellow",  "M. boulderensis"="#fecc5c", "C. pellucida"="#fd8d3c","M. sanguinipes"="#e31a1c")  
repro2$Species= factor(repro2$Species, levels=c("A. clavatus","M. boulderensis","C. pellucida","M. sanguinipes"), ordered=TRUE)

ov.plot=ggplot(data=repro2, aes(x=Elevation, y = Novarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=Novarioles-Novarioles.se, ymax=Novarioles+Novarioles.se), width=.1)+
  scale_colour_manual(values = sp.col)

clutch.plot=ggplot(data=repro2, aes(x=Elevation, y = clutch,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=clutch-clutch.se, ymax=clutch+clutch.se), width=.1)+
  scale_colour_manual(values = sp.col)

egg.mass.plot=ggplot(data=repro2, aes(x=Elevation, y = egg.mass*1000,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ xlim(1550,3530)+
  geom_errorbar(aes(ymin=egg.mass-egg.mass.se, ymax=egg.mass+egg.mass.se), width=.1)+ylab("Egg mass (mg)")+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_colour_manual(values = sp.col)

clutch.mass.plot=ggplot(data=repro2, aes(x=Elevation, y = clutch.mass,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ xlim(1550,3530)+
  geom_errorbar(aes(ymin=clutch.mass-clutch.mass.se, ymax=clutch.mass+clutch.mass.se), width=.1)+ylab("Clutch mass (g)")+xlab("Elevation (m)")+
  scale_colour_manual(values = sp.col)

Fov.plot=ggplot(data=repro2, aes(x=Elevation, y = NFovarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ xlim(1550,3530)+
  geom_errorbar(aes(ymin=NFovarioles-NFovarioles.se, ymax=NFovarioles+NFovarioles.se), width=.1)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+ylab("N functional ovarioles")+
  scale_colour_manual(values = sp.col)

PFov.plot=ggplot(data=repro2, aes(x=Elevation, y = PFovarioles,color=Species))+geom_point()+geom_line()+theme_bw()+ theme(legend.position="none")+ 
  geom_errorbar(aes(ymin=PFovarioles-PFovarioles.se, ymax=PFovarioles+PFovarioles.se), width=.1)+
  scale_colour_manual(values = sp.col)

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

#Length, CTmax, CTmin
#FOV, Egg.mass, Clutch.mass
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

#align axes
#Fov.plot, egg.mass.plot, clutch.mass.plot
gFov <- ggplotGrob(Fov.plot)
gegg.mass <- ggplotGrob(egg.mass.plot)
gclutch.mass <- ggplotGrob(clutch.mass.plot)
maxWidth = grid::unit.pmax(gFov$widths[2:5], gegg.mass$widths[2:5], gclutch.mass$widths[2:5])
gFov$widths[2:5] <- as.list(maxWidth)
gegg.mass$widths[2:5] <- as.list(maxWidth)
gclutch.mass$widths[2:5] <- as.list(maxWidth)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperGenetics/figures/")
#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperGenetics\\figures\\")
pdf("phenotype_plot.pdf", height = 10, width = 8)

grid.arrange(size.plot, ctmax.plot, ctmin.plot+ theme(legend.position = 'none'), gFov, gegg.mass, gclutch.mass, legend, hop.plot, layout_matrix = lay, widths = c(1,1), 
             heights=c(0.4,0.4,0.5,0.1,1))

dev.off()

#------------------------

phen= pbt %>% group_by(Site,Spec) %>% summarise(PBT.m=mean(PBT, na.rm=TRUE), Ctmin.m=mean(Ctmin, na.rm=TRUE),Ctmax.m=mean(Ctmax, na.rm=TRUE),mass.m=mean(mass_g, na.rm=TRUE)  )

#size
phen= size %>% group_by(Site,Species) %>% summarise(size=mean(Mean.Femur.Length, na.rm=TRUE)  )
phen=phen[order(phen$Species),]

#repro
phen= repro %>% group_by(Site,Species) %>% summarise(Mean.Femur.Length=mean(Mean.Femur.Length, na.rm=TRUE), Clutch.Size=mean(Clutch.Size, na.rm=TRUE),Mean.Egg.Mass=mean(Mean.Egg.Mass, na.rm=TRUE),Number.Ovarioles=mean(Number.Ovarioles, na.rm=TRUE),Number.Functioing.Ovarioles=mean(Number.Functioing.Ovarioles, na.rm=TRUE),Proportion.Functioning.Ovarioles=mean(Proportion.Functioning.Ovarioles, na.rm=TRUE),clutch.weight.g=mean(clutch.weight.g, na.rm=TRUE)  )
phen=phen[order(phen$Species, phen$Site),]

write.csv(phen, "PhenDat.csv" )


