# LOAD DATA BY RUNNING PhenFigs_20May2018.R

require(nlme)
require(lme4)
require(MuMIn)

#--------------
#STATS

#Fig 1. Elevation patterns

#GDDS~elevation
clim2= clim1[,c("Site","Year","elevation","period","Cdd_seas") ]
clim2= na.omit(clim2)
mod1= lme(Cdd_seas~elevation*period,random=~1|Year , data=clim2)

#---
#Phenology~elevation

dat1= as.data.frame( dat[,c("doy_adult", "gdd_adult", "elevation","species","period","site","year")] )
dat1= na.omit(dat1)

#ordinal
mod1= lme(doy_adult~elevation*species,random=~1|year , data=dat1)
#gdd
mod2= lme(gdd_adult~elevation*species,random=~1|year , data=dat1)

#----------------
#Fig 2. DI

di.dat= dat[,c("doy_adult", "gdd_adult","DI","site","year","elev","period","cdd_seas","cdd_sumfall","ordinal","species","Cdd_siteave") ]
di.dat= na.omit(di.dat)

#make random var
di.dat$spsiyr= paste(di.dat$species, di.dat$site, di.dat$year, sep="")
  
#ordinal
mod1= lme(DI~ordinal*species*period*cdd_seas*elev, random=~1|spsiyr, data=di.dat)
mod1= lme(DI~ordinal*species*cdd_seas*elev, random=~1|spsiyr, data=di.dat)
mod1= lme(DI~ordinal*species*period*Cdd_siteave, random=~1|spsiyr, data=di.dat)
#gdd
mod2= lme(DI~cdd_sumfall*species*period*cdd_seas*elev, random=~1|spsiyr, data=di.dat)
mod2= lme(DI~cdd_sumfall*species*period*Cdd_siteave, random=~1|spsiyr, data=di.dat)
#FIX ******

#divide by elev
di.dat2= subset(di.dat, di.dat$site=="C1")

#ordinal
mod1= lme(DI~ordinal*species*period*cdd_seas, random=~1|spsiyr, data=di.dat2)
#gdd
mod2= lme(DI~cdd_sumfall*species*period*cdd_seas, random=~1|spsiyr, data=di.dat2)

#--------
## STATS FOR PAPER
#divide by species
di.dat2= subset(di.dat, di.dat$species==specs[6] )
#make elevation continuous
di.dat2$elev= as.numeric(as.character(di.dat2$elev))

#dropped period from models
#ordinal
#mod1= lme(DI~ordinal+cdd_seas+elev+ordinal:cdd_seas+ordinal:elev+cdd_seas:elev, random=~1|spsiyr, data=di.dat2)
#mod1= lme(DI~ordinal*cdd_seas*elev, random=~1|spsiyr, data=di.dat2)
mod1= lme(DI~poly(ordinal, order=2)+ poly(ordinal, order=2):cdd_seas+poly(ordinal, order=2):elev+poly(ordinal, order=2):cdd_seas:elev, random=~1|spsiyr, data=di.dat2)

di.dat2$species[1]
anova(mod1, type="marginal")
#ms[,c("denDF","F-value","p-value")]
summary(mod1)

#visualize
di.plot= ggplot(data=di.dat2, aes(x=ordinal, y = DI, color=Cdd_siteave, group=spsiyr, linetype=factor(elev) ))+
  geom_smooth(method=loess, se=FALSE)+
  scale_colour_gradientn(colours =matlab.like(10))

#---
#gdd
#divide by species
di.dat2= subset(di.dat, di.dat$species==specs[5] )
#make elevation continuous
di.dat2$elev= as.numeric(as.character(di.dat2$elev))
mod2= lme(DI~poly(cdd_sumfall, order=3)+ poly(cdd_sumfall, order=3):cdd_seas+poly(cdd_sumfall, order=3):elev+poly(cdd_sumfall, order=3):cdd_seas:elev, random=~1|spsiyr, data=di.dat2)

di.dat2$species[1]
#summary(mod1)
anova(mod2, type="marginal")

#----
#divide by elevation and species
di.dat2= subset(di.dat, di.dat$species=="Melanoplus boulderensis" & di.dat$site=="C1")

#ordinal
mod1= lme(DI~ordinal*period*cdd_seas, random=~1|spsiyr, data=di.dat2)
#gdd
mod2= lme(DI~cdd_sumfall*period*cdd_seas, random=~1|spsiyr, data=di.dat2)

#----------------
#Fig 3. Adult Phenology

#aggregate to sp site year
dat.ave=dat
dat.ave$year= as.numeric(as.character(dat.ave$year))

dat.ave= aggregate(dat.ave, list(dat$spsiteyear, dat$species, dat$site, dat$period),FUN=mean, na.rm=TRUE)
dat.ave$species= dat.ave$Group.2
dat.ave$site= dat.ave$Group.3
dat.ave$period= dat.ave$Group.4

dat.s= dat.ave[,c("site","year","elevation","period","cdd_seas","doy_adult","gdd_adult","species") ]
dat.s= na.omit(dat.s)

#ordinal
mod1= lm(doy_adult~cdd_seas*species*elevation*period , data=dat.s)
#gdd
mod2= lm(gdd_adult~cdd_seas*species*elevation*period, data=dat.s)

#---
#divide by species
dat.ss= subset(dat.s, dat.s$species==specs[5])

#ordinal
mod1= lm(doy_adult~cdd_seas*elevation*period , data=dat.ss)
mod1= lm(doy_adult~cdd_seas*elevation, data=dat.ss)
mod1= lm(doy_adult~cdd_seas*site, data=dat.ss)
#gdd
mod2= lm(gdd_adult~cdd_seas*elevation*period, data=dat.ss)
mod2= lm(gdd_adult~cdd_seas*elevation, data=dat.ss)

dat.ss$species[1]
summary(mod1)

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
