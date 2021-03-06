# LOAD DATA BY RUNNING PhenFigs_20May2018.R

require(nlme)
require(lme4)
require(MuMIn)

#--------------
#STATS

#Fig 1. Elevation patterns

#subset to years with data for all
clim2= subset(clim1, clim$Year %in% c(1958:1969, 2006:2011,2014) )

#GDDS~elevation
clim3= clim2[,c("Site","Year","elevation","period","Cdd_seas") ]
clim3= na.omit(clim3)
mod1= lme(Cdd_seas~elevation*period,random=~1|Year , data=clim3)
anova(mod1)

#---
#Phenology~elevation

dat1= as.data.frame( dat[,c("doy_adult", "gdd_adult", "elevation","species","period","site","year")] )
dat1= na.omit(dat1)
#abbreviate species names
dat1$spec= substr(dat1$species, 11,14)

#restrict to years with data across sites, 2007:2011
dat1= subset(dat1, dat1$year %in% 2007:2011)

#ordinal
mod1= lme(doy_adult~elevation*species,random=~1|year , data=dat1)
#gdd
mod2= lme(gdd_adult~elevation*spec,random=~1|year , data=dat1)

#---
dat.ave.a= aggregate(dat1, list(dat1$elevation, dat1$species, dat1$year),FUN=mean, na.rm=TRUE )
dat.ave.a$species= dat.ave.a$Group.2
dat.ave.a$year= dat.ave.a$Group.3
dat.ave.a= dat.ave.a[,c(4:7,10)]

mod1= lme(doy_adult~elevation*species,random=~1|year , data=dat.ave.a)
mod2= lme(gdd_adult~elevation*species,random=~1|year , data=dat.ave.a)

#----------------
#Fig 2. DI

di.dat= dat[,c("doy_adult", "gdd_adult","DI","site","year","elev","period","cdd_seas","cdd_sumfall","ordinal","species","Cdd_siteave") ]
di.dat= na.omit(di.dat)

#make random var
di.dat$spsiyr= paste(di.dat$species, di.dat$site, di.dat$year, sep="")
  
## STATS FOR PAPER
#TABLE 1
#divide by species
di.dat2= subset(di.dat, di.dat$species==specs[5] )
#make elevation continuous
di.dat2$elev= as.numeric(as.character(di.dat2$elev))

#ordinal
#mod1= lme(DI~poly(ordinal, order=3)+ poly(ordinal, order=3):cdd_seas+poly(ordinal, order=3):elev+poly(ordinal, order=3):cdd_seas:elev, random=~1|spsiyr, data=di.dat2)
#with main effect of period
mod1= lme(DI~poly(ordinal, order=3)+ poly(ordinal, order=3):period +poly(ordinal, order=3):cdd_seas+poly(ordinal, order=3):elev+poly(ordinal, order=3):cdd_seas:elev, random=~1|spsiyr, data=di.dat2)

di.dat2$species[1]
anova(mod1, type="marginal")
#ms[,c("denDF","F-value","p-value")]
#summary(mod1)

#visualize
di.plot= ggplot(data=di.dat2, aes(x=ordinal, y = DI, color=Cdd_siteave, group=spsiyr, linetype=factor(elev) ))+
  geom_smooth(method=loess, se=FALSE)+
  scale_colour_gradientn(colours =matlab.like(10))

#---
#gdd
#divide by species
di.dat2= subset(di.dat, di.dat$species==specs[1] )
#make elevation continuous
di.dat2$elev= as.numeric(as.character(di.dat2$elev))
mod2= lme(DI~poly(cdd_sumfall, order=3)+ poly(cdd_sumfall, order=3):period + poly(cdd_sumfall, order=3):cdd_seas+poly(cdd_sumfall, order=3):elev+poly(cdd_sumfall, order=3):cdd_seas:elev, random=~1|spsiyr, data=di.dat2)

di.dat2$species[1]
#summary(mod1)
anova(mod2, type="marginal")

#analyze by GDD and ordinal date
mod3= lme(DI~poly(cdd_sumfall, order=3)+ poly(ordinal, order=3) + poly(cdd_sumfall, order=3):period + poly(cdd_sumfall, order=3):cdd_seas+poly(cdd_sumfall, order=3):elev+poly(cdd_sumfall, order=3):cdd_seas:elev, random=~1|spsiyr, data=di.dat2)

di.dat2$species[1]
#summary(mod1)
anova(mod3, type="marginal")

#----------------
#Fig 3. Adult Phenology

dat.s= na.omit(dat.ssy)
dat.s$spec= substr(dat.s$species, 11,16)
dat.s$elev= as.factor(dat.s$elev)

#ordinal
mod1= lm(doy_adult~ (cdd_seas+spec+elevation+period)^2, data=dat.s)
#gdd
mod2= lm(gdd_adult~ (cdd_seas+spec+elevation+period)^2, data=dat.s)

anova(mod1)
anova(mod2)

#---
#divide by species
dat.ss= subset(dat.s, dat.s$species==specs[6])

#ordinal
mod1= lm(doy_adult~period+cdd_seas+elevation+cdd_seas:elevation, data=dat.ss)
#gdd
mod2= lm(gdd_adult~period+cdd_seas+elevation+cdd_seas:elevation, data=dat.ss)

dat.ss$species[1]
anova(mod1)
anova(mod2)
summary(mod1)
summary(mod2)

#---
#divide by elev
dat.ss= subset(dat.s, dat.s$elev==elevs[3])

#ordinal
#mod1= lm(doy_adult~period+cdd_seas+spec+cdd_seas:spec, data=dat.ss)
mod1= lm(doy_adult~cdd_seas+spec+cdd_seas:spec, data=dat.ss)

#gdd
#mod2= lm(gdd_adult~period+cdd_seas+spec+cdd_seas:spec, data=dat.ss)
mod2= lm(gdd_adult~cdd_seas+spec+cdd_seas:spec, data=dat.ss)

dat.ss$elev[1]
anova(mod1)
anova(mod2)
summary(mod1)
summary(mod2)

#----------------
#Fig 4. Composition plot
#No stats?

#============================
##ordinal date
mod1= lm(ordinal~cdd_yr*species+elevation+period, data=hop4)

#by species
i=1
hop.sp= subset(hop4, hop4$species=specs[i])
mod1= lm(ordinal~cdd_yr*elevation+period, data=hop.sp)
summary(mod1)

##GDD
mod1= lm(cdd_sum~cdd_yr*species+elevation+period, data=hop4)

#by species
i=1
hop.sp= subset(hop4, hop4$species=specs[i])
mod1= lm(cdd_sum~cdd_yr*elevation+period, data=hop.sp)
summary(mod1)

#-------------
#STATS
hopq= subset(hop4, hop4$quantile==50)
mod1= lm(ordinal~cdd_yr*elevation+species+period,  data=hopq)

#split by species
i=1
hopq.sp= subset(hopq, hopq$species== specs[i])
mod1= lm(ordinal~cdd_yr*elevation+period,  data=hopq.sp)
mod1= lm(cdd_sum~cdd_yr*elevation+period,  data=hopq.sp)
summary(mod1)

#------------
